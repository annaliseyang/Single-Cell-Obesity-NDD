import scanpy as sc

in_path = "/home/anna_y/data/obesity_new/rna2.AD427_ADMR.QC.3343094.batch_AD427_snRNA.h5ad"
col="Class"

adata = sc.read_h5ad(in_path)

for celltype in adata.obs[col].unique():
    name = celltype.replace(" ", "_").replace("/", "_")
    print(f"\n{celltype=}, {name=}")
    subset = adata[adata.obs[col] == celltype]
    print(subset)

    out_path = f"/home/anna_y/data/obesity_new/{col}/{name}.h5ad"
    sc.write(out_path, subset)
    print(f"Subset saved to: {out_path}")
