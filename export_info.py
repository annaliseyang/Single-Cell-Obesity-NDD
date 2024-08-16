import scanpy as sc
import sys
import os
import pandas as pd

# Load data
in_path = sys.argv[1] # .h5ad file path
adata = sc.read_h5ad(in_path)
print(adata)

# parse file name
in_dir = '/'.join(in_path.split('/')[:-1])
name = in_path.split('/')[-1].split('.')[0]
celltype = name
celltype = name.split('_Aug6_2024')[0]
print(f"Processing sample: {name}", flush=True)
print(f"Cell type: {celltype}", flush=True)
# adata.obs['sample'] = sample_names

# set up output directory if it doesn't exist
out_dir = in_dir + '/' + celltype + '/'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
print(f"Output directory: {out_dir}", flush=True)

# adata.to_df().to_csv(out_dir + f'{celltype}_expression_matrix.csv')

def export_to_csv(adata, celltype, out_dir, batch_size=10000, max_file_length=None):
    # Export the expression matrix and metadata to CSV
    adata.obs.to_csv(out_dir + f'{celltype}_obs.csv')
    adata.var.to_csv(out_dir + f'{celltype}_var.csv')

    num_cells = adata.shape[0]
    print("Number of cells:", num_cells, flush=True)

    for i in range(0, num_cells, batch_size):
        if max_file_length and max_file_length < num_cells:
            n = i // max_file_length
            filename = f'{celltype}_expression_matrix_{n}.csv'
        else:
            filename = f'{celltype}_expression_matrix.csv'

        df = adata[i:i+batch_size, :].to_df()
        print(f"i={i}, Counts exported to dataframe: {df.shape}, saving to {filename}...", flush=True)
        # print(f"i={i}, n={n}, Counts exported to dataframe: {df.shape}, saving to csv file {n}...", flush=True)

        is_first_line = i % max_file_length == 0 if max_file_length else i == 0
        df.to_csv(out_dir + filename, mode = 'w' if is_first_line else 'a', header = is_first_line, index=False)
        print(f"All counts exported to csv!", flush=True)

if __name__ == "__main__":
    export_to_csv(adata, celltype, out_dir, batch_size=10000, max_file_length=1000000)
