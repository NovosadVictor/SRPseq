import cv2 as cv
import networkx as nx
from matplotlib import pyplot as plt
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
from dna_features_viewer import GraphicFeature, GraphicRecord

from src.consts import exon_len, intron_len


def plot_gene_isoforms(gene_data, output_dir):
    gene = gene_data['gene_name']
    seq_len = exon_len * len(gene_data['exon_numbers']) + intron_len * (len(gene_data['exon_numbers']) - 1)

    imgs = []
    for j, transcript in enumerate(sorted(gene_data['transcripts'], key=lambda x: x['exon_numbers'], reverse=True)):
        exon_records = [
            GraphicFeature(
                start=(exon_len + intron_len) * max(0, exon - 1),
                end=exon_len * exon + intron_len * (exon - 1),
                strand=0,
                color="#eeeeee" if exon not in transcript['exon_numbers'] else (
                    "#95d4f3" if exon in gene_data['variable_exon_numbers'] else "#c7f464"),
                label='Exon ' + str(exon),
            )
            for i, exon in enumerate(gene_data['exon_numbers'])
        ]
        #
        record = GraphicRecord(sequence_length=seq_len, features=exon_records, feature_level_height=0)
        fig, ax = plt.subplots(1, 1, figsize=(20, 2))
        record.plot(ax=ax, figure_width=40)
        ax.set_xticklabels([])
        ax.set_title(f'Isoform {j + 1}', loc='left')
        plt.savefig(f'{output_dir}/{gene}_isoform_{j}.png', bbox_inches='tight', dpi=500)
        plt.close()
        imgs.append(cv.imread(f'{output_dir}/{gene}_isoform_{j}.png'))

    imgs = cv.vconcat(imgs)
    cv.imwrite(f'{output_dir}/{gene}_isoforms.png', imgs)


def plot_isoforms_tree(tree, output_dir):
    G = nx.DiGraph()

    node_labels = {}
    edge_labels = {}

    parents = [tree]
    lvl = 0
    while len(parents):
        cur_parents = []
        for node_id, parent in enumerate(parents):
            parent_transcripts = parent.kwargs
            node_index = f'{lvl}_{node_id}'
            parent.node_id = node_index
            G.add_node(node_index)
            node_labels[node_index] = '\n'.join([t['transcript_id'] for t in parent_transcripts])

            if lvl > 0:
                G.add_edge(parent.parent.node_id, node_index)
                edge_labels[(parent.parent.node_id, node_index)] = f'{parent.divider_exon}'

            if parent.left_child is not None:
                cur_parents += [parent.left_child, parent.right_child]

        parents = cur_parents
        lvl += 1

    plt.figure(figsize=(12, 12))
    plt.title('Isoforms tree')
    pos = graphviz_layout(G, prog='dot')
    nx.draw(G, pos, labels=node_labels, with_labels=True, arrows=True, node_size=20000)
    nx.draw_networkx_nodes(G, pos)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/isoforms_tree.png', dpi=300)
