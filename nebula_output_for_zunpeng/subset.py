import scanpy as sc
import pandas as pd

tsv_path = "/net/bmc-lab4/data/kellis/group/zunpeng/HumanBrainObsity/snRNA/Nov.2/rna3.AD427_MR_Multiome.Samples.Oct22_2024.obs.tsv"
df = pd.read_csv(tsv_path, sep='\t')
print("\nAll:", df)

# Subset by "batch"=="AD427_snRNA"
subset = df[df["batch"]=="AD427_snRNA"]
print("\nSubset:", subset)
samples = subset["Sample"].tolist()

# Subset cells by "Sample" in the subset
h5ad_path = "/net/bmc-lab4/data/kellis/group/zunpeng/HumanBrainObsity/snRNA/rna2.AD427_ADMR.QC.3343094.Jun24_2024.h5ad"
out_path = "/home/anna_y/data/obesity_new/rna2.AD427_ADMR.QC.3343094.batch_AD427_snRNA.h5ad"
adata = sc.read_h5ad(h5ad_path)
print("\nadata:", adata)
adata = adata[adata.obs["Sample"].isin(samples)]
print("\nsubsetted:", adata)

sc.write(out_path, adata)
print("Subset saved to:", out_path)
