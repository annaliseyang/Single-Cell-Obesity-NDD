import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from obesity import umap_top_n_genes_by_bmi

# bmi_groups = ['bmi_<20', 'bmi_20-25', 'bmi_25-30', 'bmi_30+']

def umap_gene(adata, gene, x_groupby='bmi_groups_0', y_groupby='Subtype', save=None):
    adata = adata.copy()
    adata = adata[adata.obs[y_groupby].notna(), :] # filter cells
    sc.pp.normalize_total(adata, target_sum=1e6)
    # sc.pp.log1p(adata)

    print("Computing UMAP...", flush=True)
    sc.tl.umap(adata, random_state=42)

    # print(f"Plotting UMAP for gene {gene} in celltype {celltype}...", flush=True)
    # sc.pl.umap(adata, color=gene, title=f"{gene} in {celltype}", save=save)
    x_groups = adata.obs[x_groupby].unique()
    y_groups = adata.obs[y_groupby].unique()
    print(f"{len(x_groups)} {x_groupby} groups: {x_groups}")
    print(f"{len(y_groups)} {y_groupby} groups: {y_groups}")

    width, height = 4, 3.5 # dimensions of each subplot
    fig, axs = plt.subplots(nrows=len(y_groups), ncols=len(x_groups), figsize=(len(x_groups) * width, len(y_groups) * height), sharex=True, sharey=True)

    for row, subtype in enumerate(y_groups):
        subset = adata[adata.obs[y_groupby] == subtype]
        x_subsets = {group: subset[subset.obs[x_groupby] == group] for group in x_groups}

        values = np.concatenate([subset[:, gene].X.toarray().flatten() for subset in x_subsets.values()])
        vmax = np.percentile(values, 99)
        if vmax == 0:
            vmax = max(values)
        vmin = 0

        print(f"Plotting UMAP for gene {gene} in {subtype} subset...", flush=True)
        for col, group in enumerate(sorted(x_groups)):
            sc.pl.umap(x_subsets[group], color=gene, vmax=vmax, vmin=vmin, ax=axs[row, col])
            axs[row, col].set_title(f"{gene}, {group}\n{subtype}")

    plt.tight_layout()
    plt.savefig(save)


if __name__ == "__main__":
    in_path = sys.argv[1] # e.g. /home/anna_y/data/write/Subclass/Ast/Ast.h5ad
    name = os.path.basename(in_path).split('.')[0] # e.g. Ast

    gene = sys.argv[2] # e.g. TBC1D3D
    y_groupby = sys.argv[3] # e.g. Subtype

    print(f'Loading Anndata file: {in_path}', flush=True)
    print(f'Name: {name}', flush=True)

    adata = sc.read_h5ad(in_path)
    adata = adata[adata.obs['bmi_lv'].notna(), :] # remove cells with missing bmi values

    # deg_results_file = sys.argv[2] # e.g. /home/anna_y/results/deg_bmi_normalized/Subclass/Ast/Ast.bmi_normalized.Clean.tsv
    deg_results_file = in_path.replace('data/write', 'results/deg_bmi_normalized').replace('.h5ad', '.bmi_normalized.Ranked.Filtered.tsv')
    print(f'Loading DEG results file: {deg_results_file}', flush=True)
    # n_top = 3 # number of top positive and negative DEGs to plot

    out_path = f'figures/deg_bmi_normalized/genes/umap_{gene}_{name}_{y_groupby}_bmi_groups.png'
    umap_gene(adata, gene, x_groupby='bmi_groups_0', y_groupby=y_groupby, save=out_path)

    out_path = f'figures/deg_bmi_normalized/genes/umap_{gene}_{name}_{y_groupby}_AD_states.png'
    umap_gene(adata, gene, x_groupby='AD_states', y_groupby=y_groupby, save=out_path)
