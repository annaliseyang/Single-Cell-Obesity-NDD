import scanpy as sc
import pandas as pd
import numpy as np

in_path = "/home/anna_y/data/write/"
out_path = "/home/anna_y/data/write/celltypes/"
filename = "AD427_ADMR_meta_Jul22_2024.h5ad"

# Load the h5ad file
adata = sc.read_h5ad(in_path + filename)

adata.obs['Class'] = adata.obs['RNA.Class.Jun21_2024']
adata.obs['Subclass'] = adata.obs['RNA.Subclass.Jun21_2024']
adata.obs['Subtype'] = adata.obs['RNA.Subtype.Jun21_2024']

adata.write(in_path + 'AD427_ADMR_Aug6_2024.h5ad')

cell_types = adata.obs['Class'].unique()
print("Cell types:", cell_types, flush=True)

for cell_type in cell_types:
    print(f"Processing {cell_type}...")
    subset = adata[adata.obs['Class'] == cell_type]
    print(f"Number of cells: {len(subset)}")
    subset.write(out_path + f"{cell_type}_Aug6_2024.h5ad")
