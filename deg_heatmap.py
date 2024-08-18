import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def normalize_bmi(adata_path):
    adata = sc.read_h5ad(adata_path)
    # normalize bmi
    bmi_normalized = (adata.obs['bmi_lv'] - adata.obs['bmi_lv'].mean()) / adata.obs['bmi_lv'].std()
    adata.obs['bmi_normalized'] = bmi_normalized
    print(f'BMI normalized:', bmi_normalized[:20])
    print(adata)
    adata.write_h5ad(adata_path)

def heatmap(adata, deg_results_file, name, n_top=50, save=None):
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)

    # Read the DEG results file
    results = pd.read_csv(deg_results_file, sep='\t')

    # Extract gene names from the results
    degs = results['gene'].tolist()[:n_top]
    print(f'Top {n_top} DEGs: {degs}')

    save = f'_{name}.png'
    sc.pl.heatmap(adata, var_names=degs, log=False, standard_scale='var', vmax=1, groupby='bmi_groups', save=save)

def heatmap_with_title(adata, deg_results_file, name, n_top=50, save=None):
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)

    # Read the DEG results file
    results = pd.read_csv(deg_results_file, sep='\t')

    # Extract gene names from the results
    degs = results['gene'].tolist()[:n_top]
    print(f'Top {n_top} DEGs: {degs}')

    save = f'heatmap_{name}.png'

    # Plot the heatmap
    heatmap_dict = sc.pl.heatmap(adata, var_names=degs, log=False, standard_scale='var',
                                vmax=1, groupby='bmi_groups', save=save, show=False)

    # Set the title
    heatmap_dict['heatmap_ax'].figure.suptitle(name)

    if save:
        heatmap_dict['heatmap_ax'].figure.savefig(save)

if __name__ == "__main__":
    in_path = sys.argv[1] # e.g. /home/anna_y/data/write/Subclass/Ast/Ast.h5ad
    # in_path = "/home/anna_y/data/write/all/AD427MR_50k/AD427MR_50k.h5ad"
    name = os.path.basename(in_path).split('.')[0] # e.g. Ast
    # normalize_bmi(in_path)

    print(f'Loading Anndata file: {in_path}')
    print(f'Name: {name}')

    adata = sc.read_h5ad(in_path)
    adata = adata[adata.obs['bmi_lv'].notna(), :]
    # deg_results_file = sys.argv[2] # e.g. /home/anna_y/results/deg_bmi_normalized/Subclass/Ast/Ast.bmi_normalized.Clean.tsv
    deg_results_file = in_path.replace('data/write', 'results/deg_bmi_normalized').replace('.h5ad', '.bmi_normalized.Ranked.Filtered.tsv')
    print(f'Loading DEG results file: {deg_results_file}')
    out_path = f'_{name}.png'
    heatmap(adata, deg_results_file, name, save=out_path)
    print(f'Heatmap saved as {out_path}')
