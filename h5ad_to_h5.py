import anndata as ad

path = "/home/anna_y/data/write/"
filename = "AD427_ADMR_Aug6_2024.h5ad"

# Load the h5ad file
adata = ad.read_h5ad(path + filename)
print(adata, flush=True)

# Save the AnnData object to an HDF5 file
adata.write(path + 'AD427_ADMR_Aug6_2024.h5')
