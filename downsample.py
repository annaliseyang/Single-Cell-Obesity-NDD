import scanpy as sc
import numpy as np

path = "/home/anna_y/data/write/"
filename = "AD427_ADMR_Aug6_2024.h5ad"

# Load the h5ad file
adata = sc.read_h5ad(path + filename)

# Downsample the dataset to 1/100th of the cells
num_cells = adata.shape[0]
sample_size = num_cells // 100

# Randomly sample 1/100th of the cells
sampled_indices = np.random.choice(num_cells, sample_size, replace=False)

# Subset the AnnData object
downsampled_adata = adata[sampled_indices, :]

# Save the downsampled dataset to a new h5ad file
downsampled_adata.write(path + 'downsampled_' + filename)
