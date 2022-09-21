import pickle
import json

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split

from src.helpers.pipeline import map_motifs_to_exons, make_exons_sf_df
from src.helpers.plots import plot_isoforms_tree, plot_gene_isoforms
from src.lr import elastic_net
from src.utils.common import predict, get_scores, add_freq_to_df, make_sure_dir_exists, aggregated_score


class Pipeline:
    def __init__(self, config, gene_data, rbp_df, isoforms_df, rbps):
        self.config = config
        self.gene_data = gene_data
        self.rbp_df = rbp_df
        self.isoforms_df = isoforms_df
        self.rbps = rbps
        self.tree = None
        self.tissue_specific = self.config.get('tissue_specific', False) and 'Tissue' in self.rbp_df
        self.tissues = []
        if self.tissue_specific:
            self.tissues = set(self.rbp_df['Tissue'])
            self.rbp_df = add_freq_to_df(self.rbp_df)

        if 'Dataset.Type' in self.rbp_df.columns:
            self.train_index = self.rbp_df[self.rbp_df['Dataset.Type'] == 'Training'].index
            self.val_index = self.rbp_df[self.rbp_df['Dataset.Type'] == 'Validation'].index
        else:
            stratify = self.rbp_df['Tissue'] if self.tissue_specific else None
            self.train_index, self.val_index = train_test_split(self.rbp_df.index, test_size=.25, stratify=stratify)

    def run(self):
        exons_motifs = map_motifs_to_exons(self.gene_data, self.rbps)
        tree = make_exons_sf_df(
            self.gene_data,
            self.rbp_df, self.isoforms_df,
            gene_exon_motifs=exons_motifs,
        )

        nodes = [tree.left_child, tree.right_child]
        while len(nodes):
            for node in nodes[::2]:
                df = node.df.loc[self.train_index]
                if len(df.columns) > 2:
                    res = elastic_net(df)
                    node.res = res
                    if self.tissue_specific:
                        node.tissue_res = {}
                        for tissue in self.tissues:
                            print(tissue)
                            tissue_df = df[df['Tissue'] == tissue]
                            tissue_df['Freq'] = 1
                            res = elastic_net(tissue_df)
                            node.tissue_res[tissue] = res
            cur_nodes = []
            for node in nodes:
                if node.left_child is not None:
                    cur_nodes += [node.left_child, node.right_child]
            nodes = cur_nodes

        with open(f'{self.config["output_dir"]}/pipeline.wb', 'wb') as res_file:
            pickle.dump(self, res_file)

        plot_gene_isoforms(self.gene_data, output_dir=self.config['output_dir'])

        self.tree = tree
        self.predict()
        self.accuracy()
        self.plot()
        self.save_res()
        plot_isoforms_tree(tree, output_dir=self.config['output_dir'])

    @staticmethod
    def load_from_file(path_to_file):
        with open(path_to_file) as class_file:
            return pickle.load(class_file)

    def predict(self):
        if self.tree is None:
            return

        nodes = [self.tree.left_child, self.tree.right_child]
        while len(nodes):
            cur_nodes = []
            for node_id, node in enumerate(nodes):
                if node_id % 2 == 0:
                    node.res['predictions.train'] = predict(node.df.loc[self.train_index], node.res['coefs'], logit=False)
                    node.res['predictions.validation'] = predict(node.df.loc[self.val_index], node.res['coefs'], logit=False)

                    for tissue in self.tissues:
                        node.tissue_res[tissue]['predictions.train'] = predict(
                            node.df.loc[self.train_index].query(f'Tissue == "{tissue}"'),
                            node.tissue_res[tissue]['coefs'],
                            logit=False,
                        )
                        node.tissue_res[tissue]['predictions.validation'] = predict(
                            node.df.loc[self.val_index].query(f'Tissue == "{tissue}"'),
                            node.tissue_res[tissue]['coefs'],
                            logit=False,
                        )
                else:
                    node.res['predictions.train'] = 1 - nodes[node_id - 1].res['predictions.train']
                    node.res['predictions.validation'] = 1 - nodes[node_id - 1].res['predictions.validation']

                    for tissue in self.tissues:
                        node.tissue_res[tissue] = {
                            'predictions.train': 1 - nodes[node_id - 1].tissue_res[tissue]['predictions.train'],
                            'predictions.validation': 1 - nodes[node_id - 1].tissue_res[tissue]['predictions.validation'],
                        }

                node.res['predictions.train.accumulative'] = node.parent.res.get('predictions.train.accumulative', 1) * node.res['predictions.train']
                node.res['predictions.validation.accumulative'] = node.parent.res.get('predictions.validation.accumulative', 1) * node.res['predictions.validation']

                for tissue in self.tissues:
                    node.tissue_res[tissue]['predictions.train.accumulative'] = node.parent.tissue_res.get(tissue, {}).get('predictions.train.accumulative', 1) * node.tissue_res[tissue]['predictions.train']
                    node.tissue_res[tissue]['predictions.validation.accumulative'] = node.parent.tissue_res.get(tissue, {}).get('predictions.validation.accumulative', 1) * node.tissue_res[tissue]['predictions.validation']

                if node.left_child is not None:
                    cur_nodes += [node.left_child, node.right_child]

            nodes = cur_nodes

    def accuracy(self):
        if self.tree is None:
            return

        nodes = [self.tree.left_child, self.tree.right_child]
        while len(nodes):
            cur_nodes = []
            for node_id, node in enumerate(nodes):
                print(node.res['predictions.train'], node.df.loc[self.train_index]['fraq'])
                train_acc = get_scores(node.res['predictions.train'], node.df.loc[self.train_index]['fraq'])
                val_acc = get_scores(node.res['predictions.validation'], node.df.loc[self.val_index]['fraq'])
                node.res['accuracy'] = {
                    'train': train_acc,
                    'validation': val_acc,
                }
                for tissue in self.tissues:
                    train_acc = get_scores(node.tissue_res[tissue]['predictions.train'], node.df.loc[self.train_index].query(f'Tissue == "{tissue}"')['fraq'])
                    val_acc = get_scores(node.tissue_res[tissue]['predictions.validation'], node.df.loc[self.val_index].query(f'Tissue == "{tissue}"')['fraq'])
                    node.tissue_res[tissue]['accuracy'] = {
                        'train': train_acc,
                        'validation': val_acc,
                    }
                if node.left_child is not None:
                    cur_nodes += [node.left_child, node.right_child]
                else: # leaf node (isoform)
                    node.res['accuracy']['train.accumulative'] = get_scores(
                        node.res['predictions.train.accumulative'],
                        node.df.loc[self.train_index]['fraq'],
                    )
                    node.res['accuracy']['validation.accumulative'] = get_scores(
                        node.res['predictions.validation.accumulative'],
                        node.df.loc[self.val_index]['fraq'],
                    )

                    for tissue in self.tissues:
                        node.tissue_res[tissue]['accuracy']['train.accumulative'] = get_scores(
                            node.tissue_res[tissue]['predictions.train.accumulative'],
                            node.df.loc[self.train_index].query(f'Tissue == "{tissue}"')['fraq'],
                        )
                        node.tissue_res[tissue]['accuracy']['validation.accumulative'] = get_scores(
                            node.tissue_res[tissue]['predictions.validation.accumulative'],
                            node.df.loc[self.val_index].query(f'Tissue == "{tissue}"')['fraq'],
                        )

            nodes = cur_nodes

    def plot(self):
        if self.tree is None:
            return

        nodes = [self.tree.left_child, self.tree.right_child]
        transcript_accuracies = {}
        while len(nodes):
            cur_nodes = []
            for node_id, node in enumerate(nodes):
                if node.left_child is not None:
                    cur_nodes += [node.left_child, node.right_child]
                else:
                    transcript = node.kwargs[0]['transcript_id']
                    transcript_accuracies[transcript] = {
                        'var.train': self.isoforms_df.loc[self.train_index][transcript],
                        'var.validation': self.isoforms_df.loc[self.val_index][transcript],
                        'train': node.res['accuracy']['train.accumulative'],
                        'validation': node.res['accuracy']['validation.accumulative'],
                        'tissue': {
                            tissue: {
                                'train': node.tissue_res[tissue]['accuracy']['train.accumulative'],
                                'validation': node.tissue_res[tissue]['accuracy']['validation.accumulative'],
                                'var.train': self.isoforms_df.loc[set(self.train_index)&set(self.rbp_df[self.rbp_df['Tissue'] == tissue].index)][transcript],
                                'var.validation': self.isoforms_df.loc[set(self.val_index)&set(self.rbp_df[self.rbp_df['Tissue'] == tissue].index)][transcript],
                            } for tissue in node.tissue_res
                        }
                    }
            nodes = cur_nodes

        make_sure_dir_exists(f"{self.config['output_dir']}/scores/")
        for score in ['cor', 'mds']:
            plt.figure(figsize=(8, 6))
            bx = sns.boxplot(
                data=self.isoforms_df[list(transcript_accuracies.keys())].div(self.isoforms_df[list(transcript_accuracies.keys())].sum(axis=1), axis=0),
            )
            ax2 = bx.twinx()
            sns.scatterplot(
                x=transcript_accuracies.keys(),
                y=[transcript_accuracies[iso]['validation'][score] for iso in transcript_accuracies],
                color='blue', s=200, ax=ax2,
            )
            sns.scatterplot(
                x=transcript_accuracies.keys(),
                y=[transcript_accuracies[iso]['train'][score] for iso in transcript_accuracies],
                color='red', s=200, ax=ax2,
            )
            ax2.legend(['Validation', 'Training'])
            bx.set_xticklabels(bx.get_xticklabels(), rotation=90)
            ax2.set_xticklabels(bx.get_xticklabels(), rotation=90)
            ax2.set_ylim(0, 1)
            bx.set_ylim(0, 1)
            plt.grid(visible=True)
            plt.tight_layout()
            plt.show()
            plt.savefig(f"{self.config['output_dir']}/scores/{score}.png", dpi=300)

            if self.tissue_specific:
                tissue_scores = [aggregated_score(transcript_accuracies, tissue) for tissue in self.tissues]
                plt.figure(figsize=(8, 6))
                ax = sns.scatterplot(
                    x=list(self.tissues),
                    y=[s['validation'][score] for s in tissue_scores],
                    color='blue', s=200,
                )
                sns.scatterplot(
                    x=list(self.tissues),
                    y=[s['train'][score] for s in tissue_scores],
                    color='red', s=200, ax=ax,
                )
                ax.legend(['Validation', 'Training'])
                plt.xticks(rotation=90)
                ax.set_ylim(0, 1)
                plt.grid(visible=True)
                plt.tight_layout()
                plt.show()
                plt.savefig(f"{self.config['output_dir']}/scores/by_tissue.{score}.png", dpi=300)

    @staticmethod
    def save_res_(node, path):
        if node is not None:
            make_sure_dir_exists(path)
            if node.res:
                if node.res.get('coefs') is not None:
                    node.res['coefs'].to_csv(f'{path}/coefs.tsv', sep='\t')
                node.df.to_csv(f'{path}/df.tsv', sep='\t')
                with open(f'{path}/transcripts.json', 'w') as f:
                    json.dump({
                        'transcripts': {
                            'parent': [t['transcript_id'] for t in node.parent.kwargs] if node.parent is not None else [],
                            'node': [t['transcript_id'] for t in node.kwargs],
                        },
                        'divider_exon': node.divider_exon,
                    }, f)
                with open(f'{path}/scores.json', 'w') as f:
                    json.dump(node.res.get('accuracy'), f)

                for tissue in getattr(node, 'tissue_res', {}):
                    make_sure_dir_exists(f'{path}/by_tissue/{tissue}/')
                    if 'coefs' in node.tissue_res[tissue]:
                        node.tissue_res[tissue]['coefs'].to_csv(f'{path}/by_tissue/{tissue}/coefs.tsv', sep='\t')
                    if 'accuracy' in node.tissue_res[tissue]:
                        with open(f'{path}/by_tissue/{tissue}/scores.json', 'w') as f:
                            json.dump(node.tissue_res[tissue]['accuracy'], f)

            Pipeline.save_res_(node.left_child, path=f'{path}/left/')
            Pipeline.save_res_(node.right_child, path=f'{path}/right/')

    def save_res(self):
        if self.tree is None:
            return

        self.save_res_(self.tree, path=f"{self.config['output_dir']}/tree/")