import scanpy as sc
import os

in_dir = "/home/anna_y/data/write/"
filename = "AD427_ADMR_Aug6_2024.h5ad"

def subset_and_save(adata, condition, value, out_dir, file_suffix):
    print(f"Subsetting {condition} == {value}...", flush=True)
    print(f"Output directory: {out_dir}", flush=True)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    print(f"Processing {value}...", flush=True)
    subset = adata[adata.obs[condition] == value, :].copy()
    print(f"Number of cells: {subset.shape[0]}", flush=True)
    subset.write(os.path.join(out_dir, f"{value}_{file_suffix}.h5ad"))

def subset_celltypes(adata, celltypes, out_dir):
    for celltype in celltypes:
        subset_and_save(adata, 'Class', celltype, out_dir, "Aug6_2024")

def subset_bmi_groups(adata, bmi_groups, out_dir):
    for group in bmi_groups:
        subset_and_save(adata, 'bmi_groups', group, out_dir, "Aug6_2024")

def subset_ad_states(adata, ad_states, out_dir):
    for state in ad_states:
        subset_and_save(adata, 'AD_states', state, out_dir, "Aug6_2024")

if __name__ == "__main__":
    # Load the h5ad file
    adata = sc.read_h5ad(in_dir + filename)
    print(adata, flush=True)

    # celltypes = adata.obs['Class'].unique()
    # subset_celltypes(adata, celltypes, out_dir=os.path.join(in_dir, "celltypes"))

    # bmi_groups = ['bmi_<20', 'bmi_20-25', 'bmi_25-30', 'bmi_30+']
    bmi_groups = ['bmi_25-30', 'bmi_30+']
    subset_bmi_groups(adata, bmi_groups, out_dir=os.path.join(in_dir, "bmi_groups"))

    # ad_states = ['earlyAD', 'lateAD', 'nonAD']
    # subset_ad_states(adata, ad_states, out_dir=os.path.join(in_dir, "ad_states"))
