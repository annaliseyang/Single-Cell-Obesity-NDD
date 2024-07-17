import scanpy as sc
import os
from utils import *
from preprocessing import *
import numpy as np

init_settings()

# results_file = f"/home/anna_y/data/write/{DATASET}.h5ad"
adata = load_data(DATASET, PATH + FILENAME)
print("Data loaded!", flush=True)
print(adata, flush=True)

# plotting function
# name, func = 'TSNE', log(sc.pl.tsne)
# name, func = 'PCA', log(sc.pl.pca)
name, func = 'UMAP', log(sc.pl.umap)

variables = ['leiden_clust', 'bmi_lv', 'log1p_total_counts', 'log1p_n_genes_by_counts', 'region', 'batch']
highest_expr_genes = ['MALAT1', 'PCDH9', 'KCNIP4', 'CADM2', 'IL1RAPL1', 'DLG2', 'МТ-СО3', 'МТ-СО2', 'RBFOX1', 'SNHG14', 'CSMD1', 'CNTNAP2', 'МТ-АТР6', 'MAGI2', 'NRXN3', 'PTPRD', 'LSAMP', 'ADGRB3', 'NRG3', 'MT-CO1']
cell_types = ['RNA.Subtype.Jun21_2024', 'RNA.Class.Jun21_2024', 'RNA.Subclass.Jun21_2024']

def plot(adata, key):
    if f'{name.lower()}_{DATASET}_{key}.png' in os.listdir('/home/anna_y/scRNA/figures/'):
        print(f"{name} for {key} already exists!", flush=True)
        return

    print(f"Generating {name} for {key}...", flush=True)
    try:
        func(adata, color=key, vmax='p99', save=f'_{DATASET}_{key}.png')
        print(f"{name} saved for {key}!", flush=True)
    except Exception as e:
        print(f"Error saving {name} for {key}: {e}", flush=True)

def plot_each(adata, keys):
    for key in keys:
        if f'{name.lower()}_{DATASET}_{key}.png' in os.listdir('/home/anna_y/scRNA/figures/'):
            print(f"{name} for {key} already exists!", flush=True)
            continue

        print(f"Generating {name} for {key}...", flush=True)
        try:
            func(adata, color=key, vmax='p99', save=f'_{DATASET}_{key}.png')
            print(f"{name} saved for {key}!", flush=True)
        except Exception as e:
            print(f"Error saving {name} for {key}: {e}", flush=True)

def plot_all(adata, keys):
    keys = [key for key in keys if key in adata.var_names]
    try:
        func(adata, color=keys, vmax='p99', save=f'_{DATASET}_highest_expr_genes.png')
        print(f"{name} saved for highest_expr_genes!", flush=True)
    except Exception as e:
        print(f"Error saving {name} for highest_expr_genes: {e}", flush=True)


def plot_celltype(adata, celltype: str, key: str):
    subset = adata[adata.obs['RNA.Class.Jun21_2024'] == celltype, :]
    try:
        func(subset, color=key, save=f'_{DATASET}_{celltype}_{key}.png')
        print(f"{name} saved for {celltype}!", flush=True)
    except Exception as e:
        print(f"Error saving {name} for {celltype}: {e}", flush=True)

def plot_celltypes(adata):
    for type in adata.obs['RNA.Class.Jun21_2024'].unique():
        plot_celltype(adata, type, 'RNA.Subclass.Jun21_2024')


if __name__ == "__main__":
    # plot_all(adata, highest_expr_genes)
    plot_celltypes(adata)
