from plot import *

# if __name__ == "__main__":
# adata = load_data(f'/home/anna_y/data/write/processed_{DATASET}_Jun17_2024.h5ad')
g1 = adata[adata.obs['bmi_lv'] > 30, :]
g2 = adata[(adata.obs['bmi_lv'] > 25) & (adata.obs['bmi_lv'] < 30), :]
g3 = adata[(adata.obs['bmi_lv'] > 20) & (adata.obs['bmi_lv'] < 25), :]
g4 = adata[adata.obs['bmi_lv'] < 20, :]

def categorize_bmi(adata, num_bins=4):
    adata.obs['bmi_groups'] = pd.cut(adata.obs['bmi_lv'], bins=[0, 20, 25, 30, np.inf], labels=['bmi_<20', 'bmi_20-25', 'bmi_25-30', 'bmi_30+'])
    return adata

groups = [g1, g2, g3, g4]

highest_expr_genes = log(sc.pl.highest_expr_genes)
obesity_genes = [tup[-1] for tup in adata.uns['rank_genes_groups']['names']][:20]
anti_obesity_genes = [tup[0] for tup in adata.uns['rank_genes_groups']['names']][:20]

def highest_expr_genes_by_bmi(adata, n_top=20, save=None):
    for group, name in zip(groups, ['bmi_30+', 'bmi_25-30', 'bmi_20-25', 'bmi_<20']):
        print(f"Number of cells in group {group.obs['bmi_lv'].unique()[0]}: {len(group)}")
        # plot highest expressed genes
        highest_expr_genes(group, n_top=20, save=f'_{DATASET}_{name}.png')

        plot_celltypes(group)

        sc.pl.umap(group, color='RNA.Class.Jun21_2024', save=f'_{DATASET}_{name}_RNA.Class.Jun21_2024.png')
        # sc.pl.umap(group, color='RNA.Subtype.Jun21_2024', save=f'_{DATASET}_RNA.Subtype.Jun21_2024.png')
        # sc.pl.umap(group, color='RNA.Subclass.Jun21_2024', save=f'_{DATASET}_RNA.Subclass.Jun21_2024.png')

if __name__ == "__main__":
    adata_bmi = adata[adata.obs['bmi_lv'].notna(), :]
    sc.pl.dotplot(adata, obesity_genes, groupby='bmi_groups', standard_scale='var', save=f'{DATASET}_obesity_genes_bmi_groups.png')
    sc.pl.dotplot(adata, anti_obesity_genes, groupby='bmi_groups', standard_scale='var', save=f'{DATASET}_anti_obesity_genes_bmi_groups.png')
    sc.pl.dotplot(adata, obesity_genes, groupby='RNA.Class.Jun21_2024', standard_scale='var', save=f'{DATASET}_obesity_genes_RNA.Class.Jun21_2024.png')

    sc.pl.rank_genes_groups_heatmap(adata_bmi, n_genes=5, groupby='RNA.Class.Jun21_2024', standard_scale='var', save='_obesity_DE_RNA.Class.Jun21_2024.png')
    sc.pl.rank_genes_groups_heatmap(adata_bmi, n_genes=5, groupby='bmi_groups', standard_scale='var', save='_obesity_DE_bmi_groups.png')
