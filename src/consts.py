import os
import json


muted_columns = ['Tissue', 'fraq', 'Group', 'Dataset.Type', 'Freq']

base_dir = os.path.dirname(os.path.abspath(__file__))
genes_data = json.load(open('/huge/bulk/ENSEMBLE/genes.json', 'r'))

exon_len = 10
intron_len = 15
