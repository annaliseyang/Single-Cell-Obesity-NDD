import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from obesity import umap_top_n_genes_by_bmi

bmi_groups = ['bmi_<20', 'bmi_20-25', 'bmi_25-30', 'bmi_30+']

def deg_umap(adata, deg_results_file, n_top=3, save=None):
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e6)
    # sc.pp.log1p(adata)

    deg_results = pd.read_csv(deg_results_file, sep='\t')
    # deg_results = deg_results.sort_values(by='log2FC', ascending=False)

    pos_degs = deg_results[deg_results['log2FC'] > 0]
    neg_degs = deg_results[deg_results['log2FC'] < 0]
    print(pos_degs.head(n_top))
    print(neg_degs.head(n_top))

    degs = pos_degs['gene'].tolist()[:n_top] + neg_degs['gene'].tolist()[:n_top]
    print(f'Top {n_top} DEGs: {degs}')

    print("Computing UMAP...", flush=True)
    sc.tl.umap(adata, random_state=42)

    # subset adata by bmi groups
    bmi_subsets = {group: adata[adata.obs['bmi_groups'] == group] for group in bmi_groups}
    print(f'BMI subsets: {bmi_subsets.keys()}')

    umap_top_n_genes_by_bmi(bmi_subsets, colors=degs, n_top=len(degs), save=save)


if __name__ == "__main__":
    in_path = sys.argv[1] # e.g. /home/anna_y/data/write/Subclass/Ast/Ast.h5ad
    name = os.path.basename(in_path).split('.')[0] # e.g. Ast

    print(f'Loading Anndata file: {in_path}', flush=True)
    print(f'Name: {name}', flush=True)

    adata = sc.read_h5ad(in_path)
    adata = adata[adata.obs['bmi_lv'].notna(), :] # remove cells with missing bmi values

    # deg_results_file = sys.argv[2] # e.g. /home/anna_y/data/results/deg_bmi_normalized/Subclass/Ast/Ast.bmi_normalized.Clean.tsv
    deg_results_file = in_path.replace('data/write', 'results/deg_bmi_normalized').replace('.h5ad', '.bmi_normalized.Ranked.Filtered.tsv')
    print(f'Loading DEG results file: {deg_results_file}', flush=True)
    n_top = 3 # number of top positive and negative DEGs to plot

    out_path = f'/home/anna_y/data/results/figures/deg_bmi_normalized/celltypes/umap_{name}_bmi_groups.png'
    deg_umap(adata, deg_results_file, n_top=n_top, save=out_path)
    print(f'UMAP saved to {out_path}')
