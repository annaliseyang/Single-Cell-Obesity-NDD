import anndata as ad

path = "/home/anna_y/data/write/"
filename = "AD427_ADMR_meta_Jul22_2024.h5ad"

# Load the h5ad file
adata = ad.read_h5ad(path + filename)
print(adata, flush=True)

# Save the AnnData object to an HDF5 file
adata.write(path + 'AD427_ADMR_Aug6_2024.h5')


# # Export the expression matrix and metadata to CSV
# adata.to_df().to_csv(path + 'expression_matrix.csv')
# adata.obs.to_csv(path + 'metadata.csv')
# adata.var.to_csv(path + 'features.csv')

# print("Data exported to CSV", flush=True)
