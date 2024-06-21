import pandas as pd
import gseapy as gp
import numpy as np
import os
import json
from differential_expression.code._constants import *

METRIC = "NOM p-val"

for TERM in TERMS:

    if TERM == "BP":
        GENE_SET = 'GO-biological_process'
        top_number = 100
    elif TERM == "MF":
        GENE_SET = 'GO-molecular_function'
        top_number = 100
    elif TERM == "DB":
        GENE_SET = 'Dictybase-Phenotypes'
        top_number = 50

    data_path = '../data/combined_df_entrez.csv'
    gene_sets_path = f'../data/genesets/{GENE_SET}_mod.gmt'


    data_t = pd.read_csv(data_path)
    gene_set_descriptions = {}
    with open(gene_sets_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            gene_set_name = parts[0]
            gene_set_description = parts[1] if len(parts) > 1 else "No description"
            gene_set_descriptions[gene_set_name] = gene_set_description
    available_cores = os.cpu_count()
    for hour in hours:
        for strain in strains:
            current_data = data_t[(data_t['hour'] == hour) & ((data_t['strain'] == strain) | (data_t['strain'] == 'AX4'))]
            if current_data.empty:
                continue  
            data = current_data.T
            data = data.drop(['replicate', 'hour', 'strain'])
            metadata = current_data['strain']
            data = data.astype(float)
            gs_res = gp.gsea(data=data,
                            gene_sets=gene_sets_path,
                            cls=metadata.tolist(),
                            permutation_num=400,
                            no_plot=True,
                            method='signal_to_noise',
                            threads=available_cores,
                            format='png')
            gs_res.res2d['Description'] = gs_res.res2d['Term'].map(gene_set_descriptions)
            top_gene_sets = gs_res.res2d.sort_values(METRIC, ascending=True).head(top_number)
            top_gene_sets['Rank'] = np.arange(1, len(top_gene_sets) + 1)
            top_gene_sets = top_gene_sets[['Rank', 'Term', 'Description', METRIC]]
            top_gene_sets['Description'] = top_gene_sets['Description'].str[1:-1]
            top_gene_sets = top_gene_sets.to_dict(orient='records')
            for entry in top_gene_sets:
                entry['adj.pValue'] = entry.pop('NOM p-val')
            out_file = f'../data/pairwise/gsea/{TERM}/DE_{hour}_{strain}.json'
            os.makedirs(os.path.dirname(out_file), exist_ok=True)
            with open(out_file, 'w') as outfile:
                json.dump(top_gene_sets, outfile, indent=4)
