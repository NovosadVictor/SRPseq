import pandas as pd
import sys
import os
import re
import json
from functools import partial
from shutil import copyfile
import matplotlib.pyplot as plt
import numpy as np

from src.consts import genes_data, base_dir, muted_columns
from src.helpers.plots import plot_gene_isoforms
from src.tree import TranscriptsTreeNode
from src.utils.common import find_substring_occurrences, intersect_dfs, make_sure_dir_exists


def load_rbp_data():
    return pd.read_csv('/huge/bulk/TCGA/TCGA-COMBINED/combined/sfs_FPKM.tsv', sep='\t', index_col=0).T


def load_isoforms(gene_name):
    return 2**pd.read_csv(f'/huge/bulk/TCGA/TCGA-COMBINED/isoforms/by_gene/{gene_name}_isoform_FPKM.tsv', sep='\t', index_col=0) - 1


def load_rbps():
    return pd.read_csv(f'{base_dir}/../data/new_splicing_factors_symbols.tsv', sep='\t', index_col=0)


def get_first_variable_exon(transcripts, starting=-1):
    all_exons = list(sorted(set(e for transcript in transcripts for e in transcript['exon_numbers'] if e > starting)))
    for e in all_exons:
        for transcript in transcripts:
            if e not in transcript['exon_numbers']:
                return e

    return None


def set_variable_exons(gene_data):
    common_exons = [e for e in gene_data['exon_numbers'] if all([e in transcript['exon_numbers'] for transcript in gene_data['transcripts']])]
    gene_data['variable_exon_numbers'] = [e for e in gene_data['exon_numbers'] if e not in common_exons]
    gene_data['variable_exons'] = [e for e in gene_data['exons'] if e['exon_number'] not in common_exons]
    for transcript in gene_data['transcripts']:
        transcript['variable_exon_numbers'] = [e for e in transcript['exon_numbers'] if e not in common_exons]
        transcript['variable_exons'] = [e for e in transcript['exons'] if e['exon_number'] not in common_exons]

    return gene_data


def map_exons_to_numbers(gene_data):
    gene_data['exon_numbers'] = [int(e['exon_number']) for e in gene_data['exons']]
    for i in range(len(gene_data['transcripts'])):
        gene_data['transcripts'][i]['exon_numbers'] = [int(e['exon_number']) for e in gene_data['transcripts'][i]['exons']]

    return gene_data


def load_data(gene_name):
    return map_exons_to_numbers(genes_data.copy().get(gene_name))


def filter_columns_by_expression(df, tresh_mean, tresh_var):
    return df.drop(columns=[
        c for c in df
        if np.issubdtype(df[c].dtype, np.number) and (df[c].mean() < tresh_mean or df[c].var() < tresh_var)
    ])


def make_exon_sf_df(sfs_df, isoforms_df, gene_exon_motifs, exon_number, node_isoforms, parent_isoforms, tr_low=1.0, tr_high=0.0):
    exon_sfs = set(gene_exon_motifs[
       (gene_exon_motifs['Number'].astype(int) == exon_number)
       & (
               (gene_exon_motifs['Pos'] == 'Exon')
               | ((gene_exon_motifs['Loc.Percent'] >= tr_high) | (gene_exon_motifs['Loc.Percent'] <= tr_low))
       )
    ].index)
    #
    exon_df = sfs_df[set(exon_sfs) & set(sfs_df.columns)]
    exon_df.loc[:, 'fraq'] = isoforms_df[node_isoforms].sum(axis=1) / isoforms_df[parent_isoforms].sum(axis=1)
    exon_df.loc[:, 'fraq'] = (exon_df['fraq'] * (len(exon_df) - 1) + 0.5) / len(exon_df)
    exon_df.loc[:, set(sfs_df.columns) - set(exon_df.columns)] = sfs_df[set(sfs_df.columns) - set(exon_df.columns)]
    #
    return exon_df


def make_transcripts_tree(transcripts, exons):
    transcripts_tree = TranscriptsTreeNode(kwargs=transcripts)
    parents = [transcripts_tree]
    for lvl in range(len(exons)):
        cur_parents = []
        for parent in parents:
            parent_transcripts = parent.kwargs
            if len(parent_transcripts) > 1:
                variable_exon = get_first_variable_exon(parent_transcripts, starting=parent.divider_exon)
                if variable_exon is not None:
                    lhs_node, rhs_node = [], []
                    for transcript in parent_transcripts:
                        if variable_exon in transcript['variable_exon_numbers']:
                            lhs_node.append(transcript)
                        else:
                            rhs_node.append(transcript)
                    lhs_node = TranscriptsTreeNode(kwargs=lhs_node, parent=parent, divider_exon=variable_exon)
                    rhs_node = TranscriptsTreeNode(kwargs=rhs_node, parent=parent, divider_exon=variable_exon)
                    parent.set_children(left=lhs_node, right=rhs_node)
                    cur_parents += [lhs_node, rhs_node]
        parents = cur_parents

    return transcripts_tree


def make_exons_sf_df(gene_data, sfs_df, isoforms_df, gene_exon_motifs):
    transcripts = gene_data['transcripts']
    transcripts = [t for t in transcripts if t['transcript_id'] in isoforms_df.columns]
    exons = gene_data['variable_exon_numbers']

    transcripts_tree = make_transcripts_tree(transcripts, exons)

    nodes = [transcripts_tree.left_child, transcripts_tree.right_child]
    while len(nodes):
        for node_id, node in enumerate(nodes):
            if node_id % 2 == 0:
                parent_transcripts = node.parent.kwargs
                node_transcripts = node.kwargs
                node.df = make_exon_sf_df(
                    sfs_df, isoforms_df, gene_exon_motifs,
                    node.divider_exon,
                    [t['transcript_id'] for t in node_transcripts],
                    [t['transcript_id'] for t in parent_transcripts],
                )
            else:
                node.df = nodes[node_id - 1].df.copy()
                node.df['fraq'] = 1 - node.df['fraq']

        cur_nodes = []
        for node in nodes:
            if node.left_child is not None:
                cur_nodes += [node.left_child, node.right_child]
        nodes = cur_nodes

    return transcripts_tree


def map_motifs_to_exons(gene_data, motifs_data):
    motifs = motifs_data.copy()
    exons = pd.DataFrame(gene_data['exons']).drop_duplicates()
    variable_exons = pd.DataFrame(gene_data['variable_exons']).drop_duplicates()
    if variable_exons.iloc[0]['exon_number'] != 1:
        start = exons[exons['exon_number'] == variable_exons.iloc[0]['exon_number'] - 1].iloc[0]['end'] - 1
    else:
        start = 0
    end = variable_exons.iloc[-1]['end']
    gene_seq = gene_data['sequence'][start: end]
    variable_exons.loc[:, ['start', 'end']] = variable_exons.loc[:, ['start', 'end']] - start
    motifs['in gene'] = motifs['Motif'].apply(lambda x: x in gene_seq)
    motifs = motifs[motifs['in gene'] == True]

    motif_exons = pd.DataFrame()
    for motif in set(motifs['Motif']):
        motif_locs = find_substring_occurrences(motif, gene_seq)
        motif_chunk = []
        for loc in motif_locs:
            motif_chunk.append(find_nearest_exon(loc, variable_exons))
        motif_chunk = pd.DataFrame(motif_chunk)
        motif_chunk['Motif'] = [motif] * len(motif_chunk)
        motif_exons = pd.concat([motif_exons, motif_chunk], axis=0)

    return motif_exons.merge(motifs_data, on='Motif').set_index('Gene')


def find_nearest_exon(loc, exons):
    inside_exon = exons[(exons['start'] <= loc) & (exons['end'] >= loc)]
    if not inside_exon.empty:
        inside_exon = inside_exon.iloc[0]
        return {
            'Pos': 'Exon',
            'Number': inside_exon['exon_number'],
            'Loc.Absolute': loc,
            'Loc.Relative': loc - inside_exon['start'],
            'Loc.Percent': (loc - inside_exon['start']) / (inside_exon['end'] - inside_exon['start']),
        }

    idx = exons[(exons['start'] >= loc)]['start'].idxmin()
    nearest_exon = exons.loc[idx]
    if nearest_exon.name == 0:
        prev_end = 0
    else:
        prev_end = exons.iloc[exons.index.get_loc(nearest_exon.name) - 1]['end']

    return {
        'Pos': 'Intron',
        'Number': nearest_exon['exon_number'],
        'Loc.Absolute': loc,
        'Loc.Relative': loc - prev_end,
        'Loc.Percent': (loc - prev_end) / (nearest_exon['start'] - prev_end),
    }


def load_config_and_input_data(config_path):
    """Load configuration file and input data
    Parameters
    ----------
    config_path : string
        Path to config file (json).
    Returns
    -------
    dict, dict, pd.DataFrame, pd.DataFrame, pd.DataFrame
    """
    #
    print('Loading config...')
    try:
        config_file = open(config_path, 'r')
        config = json.load(config_file)
    except:
        print('Cannot open configuration file', file=sys.stderr)
        sys.exit(1)
    #
    # Paths are absolute or relative to config file
    config_dirname = os.path.dirname(config_path)
    #
    gene = config['gene']
    if config.get('rbp_data_path') and os.path.isfile(config['rbp_data_path']):
        rbp_df = pd.read_csv(os.path.join(config_dirname, config['rbp_data_path']), sep='\t', index_col=0)
    else:
        rbp_df = load_rbp_data()
    rbp_df = filter_columns_by_expression(
        rbp_df,
        tresh_mean=config.get('rbps_tresh_mean', 1),
        tresh_var=config.get('rbps_tresh_var', 3),
    ).rename(columns={'Cancer': 'Tissue'})
    #
    if config.get('isoforms_data_path') and os.path.isfile(config['isoforms_data_path']):
        isoforms_df = pd.read_csv(os.path.join(config_dirname, config['isoforms_data_path']), sep='\t', index_col=0)
    else:
        isoforms_df = load_isoforms(gene)
    isoforms_df = filter_columns_by_expression(
        isoforms_df,
        tresh_mean=config.get('isoforms_tresh_mean', 1),
        tresh_var=config.get('isoforms_tresh_var', 10),
    )
    if config.get('rbps_path') and os.path.isfile(config['rbps_path']):
        rbps = pd.read_csv(os.path.join(config_dirname, config['rbps_path']), sep='\t', index_col=0)
    else:
        rbps = load_rbps()
    #
    rbp_df, isoforms_df = intersect_dfs([rbp_df, isoforms_df])
    #
    gene_data = load_data(gene)
    gene_data['transcripts'] = [t for t in gene_data['transcripts'] if t['transcript_id'] in isoforms_df.columns]
    gene_data = set_variable_exons(gene_data)
    gene_data['sequence'] = gene_data['sequence'].replace('T', 'U')
    #
    config['output_dir'] = make_sure_dir_exists(os.path.join(config_dirname, config['output_dir']))
    copyfile(config_path, os.path.join(config['output_dir'], 'config.json'))
    #
    print('Loaded config...')
    #
    return config, gene_data, rbp_df.astype({col: np.float64 for col in rbp_df.columns if col not in muted_columns}), isoforms_df, rbps
