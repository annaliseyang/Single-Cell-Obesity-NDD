import scanpy as sc
import os
import sys

in_dir = "/home/anna_y/data/write/"
filename = "AD427_ADMR_Aug6_2024.h5ad" # full dataset

group_by = sys.argv[1]

def subset_and_save(adata, condition, value, out_dir, file_suffix=""):
    print(f"\nSubsetting {condition} == {value}...", flush=True)
    print(f"Output directory: {out_dir}", flush=True)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    print(f"Processing {value}...", flush=True)
    subset = adata[adata.obs[condition] == value, :].copy()
    print(f"Number of cells: {subset.shape[0]}", flush=True)
    value = str(value).replace(" ", "_").replace("/", "_")  # replace spaces and slashes with underscores for file name
    if not file_suffix:
        subset.write(os.path.join(out_dir, f"{value}.h5ad"))
    else:
        subset.write(os.path.join(out_dir, f"{value}_{file_suffix}.h5ad"))

def subset_all(adata, group_by, out_dir):
    groups = adata.obs[group_by].unique()
    for group in groups:
        subset_and_save(adata, group_by, group, out_dir)

if __name__ == "__main__":
    # Load the h5ad file
    adata = sc.read_h5ad(in_dir + filename)
    print(adata, flush=True)

    # out_dir = os.path.join(in_dir, group_by)
    subset_all(adata, group_by, out_dir=os.path.join(in_dir, group_by))
