import scanpy as sc
import numpy as np

in_dir = "/home/anna_y/data/write/"
out_dir = "/home/anna_y/data/test/"
filename = "AD427_ADMR_Aug6_2024.h5ad"

# Load the h5ad file
adata = sc.read_h5ad(in_dir + filename)

# Downsample the dataset
num_cells = adata.shape[0]
sample_size = num_cells // 1000

# Randomly select indices for downsampling
sampled_indices = np.random.choice(num_cells, sample_size, replace=False)

# Subset the AnnData object
downsampled_adata = adata[sampled_indices, :]
print("Downsampled dataset:")
print(downsampled_adata, flush=True)

# Save the downsampled dataset to a new h5ad file
out_path = out_dir + 'tiny_' + filename
downsampled_adata.write(out_path)
print("Downsampled dataset saved to:", out_path, flush=True)
