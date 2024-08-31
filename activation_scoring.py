import scanpy as sc
import sys
import os

mic_genes = ["TPT1", "DUSP1", "B2M", "TREM2", "CCL2", 'APOE', "AXL", "ITGAX", "CD9", "C1QA", "C1QC", "CTSS", "CCL3","CSF3R","CX3CR1","SLC2A5","TMEM119","CD68"]

def activation_score(adata, genes, groupby='bmi_groups', save=None):
    """
    Calculate the activation score for each cell in the provided list of genes.
    """
    # adata.obs['activation_score'] = adata.X[:, adata.var_names.isin(genes)].sum(axis=1)
    sc.tl.score_genes(adata, gene_list=genes, score_name='mic_score')
    # print(adata.obs['mic_score'].head())
    avg = adata.obs['mic_score'].mean()
    sc.pl.violin(adata, keys='mic_score', groupby=groupby, save=save)
    sc.pl.heatmap(adata, var_names=genes, groupby=groupby, save=save)
    # print("average: ", avg)
    return avg

if __name__ == "__main__":
    # in_path = "/home/anna_y/data/write/Class/Oli_50k/Oli_50k.h5ad"
    in_path = sys.argv[1] # e.g. /home/anna_y/data/write/Class/Oli/Oli.h5ad
    celltype = os.path.basename(in_path).split('.')[0]

    adata = sc.read_h5ad(in_path)
    print(f'Computing UMAP of mic_genes in {celltype}...', flush=True)
    sc.tl.umap(adata, random_state=42)
    sc.pl.umap(adata, color=mic_genes, vmax='p99', save=f'_mic_genes_{celltype}.png')
    # group = adata[adata.obs['bmi_groups'] == 'bmi_20-25', :].copy()
    # adata = activation_score(group, mic_genes)
    groupby_list=['bmi_groups', 'obesity_groups', 'AD_states']
    for groupby in groupby_list:
        print("\nGrouping by: ", groupby)
        adata = adata[adata.obs[groupby].notna(), :].copy()
        activation_score(adata, mic_genes, groupby=groupby, save=f'_mic_genes_{celltype}_{groupby}.png')

        for group in adata.obs[groupby].unique().dropna():
            print(f"Processing {group}")
            subset = adata[adata.obs[groupby] == group, :].copy()
            score = activation_score(subset, mic_genes)
            print(f"mic_genes activation for {group}: {score}")
