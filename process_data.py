import scanpy as sc
import os
from preprocessing import *
import numpy as np

if __name__ == "__main__":
    init_settings()
    results_file = f"/home/anna_y/data/write/{DATASET}_new.h5ad"

    data = load_data(DATASET, results_file)
    print("Data loaded!", flush=True)
    print(data, flush=True)

    # # compute neighbors
    # neighbors = log(sc.pp.neighbors)
    # neighbors(data)
    # print("Neighbors computed!", flush=True)

    # plot highest expressed genes
    highest_expr_genes = log(sc.pl.highest_expr_genes)
    highest_expr_genes(data, n_top=20, save=f'_{DATASET}.png')

    # add the highly expressed genes annotation to the obs
    data.obs['highest_expr_genes'] = data.var_names[data.var[:, 'n_cells'] > data.obs['n_cells'].max() * 0.01]
    # print the highly expressed genes
    print("Highest expressed genes:", data.obs['highest_expr_genes'].tolist(), flush=True)

    # make UMAP
    umap = log(sc.pl.umap)
    umap(data, color='highest_expr_genes', save=f'_{DATASET}_highest_expr_genes.png')

    print(data, flush=True)

    # save results
    data.write(f'/home/anna_y/data/write/processed_{DATASET}.h5ad')

    # # annotate celltype
    # sc.tl.annotate_cell_groups(data, key='celltype', groupby='louvain')
    # print("Celltype annotation completed!", flush=True)

    # for key, val in data.var.highly_variable.items():
    #     if val == True:
    #         print(f"Highly variable gene: {key}", flush=True)
