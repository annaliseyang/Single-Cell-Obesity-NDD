import scanpy as sc
import sys
import os

# Load data
in_path = sys.argv[1]
adata = sc.read_h5ad(in_path)

# parse file name
in_dir = '/'.join(in_path.split('/')[:-1])
name = in_path.split('/')[-1].split('.')[0]
celltype = name.split('_Aug6_2024')[0]
print(f"Processing sample: {name}", flush=True)
print(f"Cell type: {celltype}", flush=True)
# adata.obs['sample'] = sample_names

# set up output directory if it doesn't exist
out_dir = in_dir + '/' + celltype + '/'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
print(f"Output directory: {out_dir}", flush=True)

# Export the expression matrix and metadata to CSV
adata.to_df().to_csv(out_dir + f'{celltype}_expression_matrix.csv')
adata.obs.to_csv(out_dir + f'{celltype}_obs.csv')
adata.var.to_csv(out_dir + f'{celltype}_var.csv')

print("Data exported to CSV", flush=True)
