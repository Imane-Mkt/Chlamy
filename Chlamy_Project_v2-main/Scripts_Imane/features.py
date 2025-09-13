import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

def add_GO_old(data):
    gene_info = pd.read_csv('/home/imokhtatif/.vscode-server/Chlamy_Project_v2-main/CreinhardtiiCC_4532_707_v6.1.annotation_info.txt', sep='\t')
    gene_GO = gene_info[['Best-hit-clamy-name', 'GO']].dropna()
    gene_GO.rename(columns={'Best-hit-clamy-name': 'mutated_genes'}, inplace=True)
    gene_GO['GO'] = gene_GO['GO'].apply(lambda x: x.split(' '))
    gene_GO.drop_duplicates(subset='mutated_genes', inplace=True)

    data = pd.merge(data, gene_GO, on='mutated_genes', how='left')
    data['GO'] = data['GO'].apply(lambda x: [] if not isinstance(x, list) and pd.isna(x) else x)
    return(data)


def add_GO(data, annotation_file):
    def uniq_list(seq):
        return list(dict.fromkeys([x for x in seq if pd.notna(x) and str(x).strip() != ""]))
    gene_info = pd.read_csv(annotation_file, sep=",", dtype=str)
    gene_info['mutated_genes'] = gene_info['Gene Name'].str.replace(r'_4532$', '', regex=True)
    cols = ['mutated_genes', 'GO ID', 'GO Description', 'Pathway Description', 'Panther ID']#, 'KEGG ID']
    gene_GO = gene_info[cols]#.dropna(subset=['GO ID'])
    gene_GO = gene_GO.groupby('mutated_genes', as_index=False).agg({
        'GO ID': uniq_list,
        'GO Description': uniq_list,
        'Pathway Description': uniq_list,
        'Panther ID': uniq_list
        #'KEGG ID': uniq_list
    })
    out = pd.merge(data, gene_GO, on='mutated_genes', how='left')
    for col in ['GO ID', 'GO Description', 'Pathway Description', 'Panther ID']:#, 'KEGG ID']:
        out[col] = out[col].apply(lambda x: [] if not isinstance(x, list) and pd.isna(x) else x)
    out['GO'] = out['GO ID']
    return out

