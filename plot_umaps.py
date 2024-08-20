import scanpy as sc
import os
from utils import *
from preprocessing import *

init_settings()

# results_file = f"/home/anna_y/data/write/{DATASET}.h5ad"
data = load_data(DATASET, PATH + FILENAME)
print("Data loaded!", flush=True)
print(data, flush=True)

# plot umaps
umap = log(sc.pl.umap)

for key in data.obs.columns:
    if f'umap_{DATASET}_{key}.png' in os.listdir('/home/anna_y/scRNA//home/anna_y/data/results/figures/'):
        print(f"UMAP for {key} already exists!", flush=True)
        continue

    print(f"Generating UMAP for {key}...", flush=True)
    try:
        umap(data, color=key, vmax='p99', save=f'_{DATASET}_{key}_p99.png')
        print(f"UMAP saved for {key}!", flush=True)
    except Exception as e:
        print(f"Error saving UMAP for {key}: {e}", flush=True)
