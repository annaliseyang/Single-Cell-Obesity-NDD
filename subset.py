import scanpy as sc
import os
import sys

# in_dir = "/home/anna_y/data/write/"
# filename = "AD427_ADMR_Aug6_2024.h5ad" # full dataset
in_path = sys.argv[1] # path to full dataset, e.g. "/home/anna_y/data/write/AD427_ADMR_Aug6_2024.h5ad"
group_by_inputs = sys.argv[2:] # e.g. "Class"

def subset_and_save(adata, condition, value, out_dir=None, file_suffix=None, out_path=None):
    print(f"\nSubsetting {condition} == {value}...", flush=True)
    print(f"Output directory: {out_dir}", flush=True)

    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)
    print(f"Processing {value}...", flush=True)
    subset = adata[adata.obs[condition] == value, :].copy()
    print(f"Number of cells: {subset.shape[0]}", flush=True)
    value = str(value).replace(" ", "_").replace("/", "_")  # replace spaces and slashes with underscores for file name
    if out_path:
        subset.write(out_path)
    elif not file_suffix:
        subset.write(os.path.join(out_dir, f"{value}.h5ad"))
    else:
        subset.write(os.path.join(out_dir, f"{value}_{file_suffix}.h5ad"))

def subset_all(adata, group_by, out_dir):
    groups = adata.obs[group_by].unique()
    for group in groups:
        subset_and_save(adata, group_by, group, out_dir)

if __name__ == "__main__":
    # Load the h5ad file
    adata = sc.read_h5ad(in_path)
    print(adata, flush=True)

    # out_dir = os.path.join(in_dir, group_by)
    in_dir = os.path.dirname(in_path)
    # for group_by in group_by_inputs:
    #     print(f"\nProcessing {group_by}...", flush=True)
    #     subset_all(adata, group_by, out_dir=os.path.join(in_dir, group_by))

    # Subset region == PFC
    out_path = in_path.replace(".h5ad", ".PFC.h5ad")
    subset_and_save(adata, "region", "PFC", out_path=out_path)
    print("PFC subset saved to", out_path, flush=True)
