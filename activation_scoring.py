import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

mic_genes = ["TPT1", "DUSP1", "B2M", "TREM2", "CCL2", 'APOE', "AXL", "ITGAX", "CD9", "C1QA", "C1QC", "CTSS", "CCL3","CSF3R","CX3CR1","SLC2A5","TMEM119","CD68"]

def activation_score(adata, genes, celltype = None, groupby='bmi_groups', score_name="activation_score", save=None):
    """
    Calculate the activation score for each cell in the provided list of genes.
    """
    # adata.obs['activation_score'] = adata.X[:, adata.var_names.isin(genes)].sum(axis=1)
    sc.tl.score_genes(adata, gene_list=genes, score_name=score_name)
    # print(adata.obs['mic_score'].head())
    avg = adata.obs[score_name].mean()
    # sc.pp.normalize_per_cell(adata, key=score_name, layer=None, copy=False)
    # sc.tl.norm(adata, key=score_name, layer=None, copy=False)
    # sc.pl.violin(adata, keys=score_name, groupby=groupby, scale='area', stripplot=False, save=save)
    # sc.pl.scatter(adata, x='bmi_lv', y=score_name, color='Subclass', save=save)
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.boxplot(x=groupby, y=score_name, data=adata.obs, hue=groupby, palette='viridis', legend=False, ax = ax, showfliers=False)
    plt.title(f"{score_name} ({celltype})")
    plt.savefig(save)
    plt.close()
    # print("average: ", avg)
    return avg

if __name__ == "__main__":
    # in_path = "/home/anna_y/data/write/Class/Oli_50k/Oli_50k.h5ad"
    in_path = sys.argv[1] # e.g. /home/anna_y/data/write/Class/Mic_Immune/Mic_Immune.h5ad
    celltype = os.path.basename(in_path).split('.')[0]

    adata = sc.read_h5ad(in_path)
    print(f'Computing UMAP of mic_genes in {celltype}...', flush=True)
    sc.tl.umap(adata, random_state=42)
    sc.pl.umap(adata, color=mic_genes, vmax='p99', save=f'_mic_genes_{celltype}.png')
    # group = adata[adata.obs['bmi_groups'] == 'bmi_20-25', :].copy()
    # adata = activation_score(group, mic_genes)
    groupby_list=['bmi_groups', 'obesity_groups', 'AD_states']
    for groupby in groupby_list:
        print("\nGrouping by: ", groupby)
        adata = adata[adata.obs[groupby].notna(), :].copy()
        activation_score(adata, mic_genes, celltype=celltype, groupby=groupby, score_name='mic_score', save=f'_mic_genes_{celltype}_{groupby}.png')

        for group in adata.obs[groupby].unique().dropna():
            print(f"Processing {group}")
            subset = adata[adata.obs[groupby] == group, :].copy()
            score = activation_score(subset, mic_genes)
            print(f"mic_genes activation for {group}: {score}")
