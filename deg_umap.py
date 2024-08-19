import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from obesity import umap_by_bmi, umap_by_groups, umap_top_n_genes_by_bmi, umap_top_n_genes_by_groups
import sys
import os

def deg_umaps(adata, name, deg_results_file, n_top=5, save=None):
    """
    Plot UMIs for top n genes in each group.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    deg_results_file (str): The path to the file containing the DEG results.
    n_top (int): The number of top genes to plot.
    """
    # cluster cell types
    adata = adata.copy()
    # sc.tl.leiden(adata)

    # Load the DEG results and get top genes
    results = pd.read_csv(deg_results_file, sep='\t')
    degs = results['gene'].tolist()[:n_top]
    print(f'Top {n_top} DEGs: {degs}')

    # Subset by bmi groups
    bmi_groups = ['bmi_<20', 'bmi_20-25', 'bmi_25-30', 'bmi_30+']
    bmi_subsets = {group: adata[adata.obs['bmi_groups'] == group] for group in bmi_groups}

    # UMIs for top n genes by bmi groups
    umap_top_n_genes_by_bmi(bmi_subsets, degs, save=save)


if __name__ == "__main__":
    in_path = sys.argv[1] # e.g. /home/anna_y/data/write/Subclass/Ast/Ast.h5ad
    name = os.path.basename(in_path).split('.')[0] # e.g. Ast

    # Load the AnnData file
    adata = sc.read_h5ad(in_path)

    # Deg UMIs
    deg_results_file = in_path.replace('data/write', 'results/deg_bmi_normalized').replace('.h5ad', '.bmi_normalized.Ranked.tsv')
    out_path = f'figures/deg_bmi_normalized/umap_{name}_bmi_groups.png'
    deg_umaps(adata, name, deg_results_file, save=out_path)
    print(f'DEG UMIs saved to {out_path}')
