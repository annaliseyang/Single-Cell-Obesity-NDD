from plot import *
import pandas as pd

# if __name__ == "__main__":
adata = load_data(f'/home/anna_y/data/write/processed_{DATASET}_Jun17_2024.h5ad')
g1 = adata[adata.obs['bmi_lv'] > 30, :]
g2 = adata[(adata.obs['bmi_lv'] > 25) & (adata.obs['bmi_lv'] < 30), :]
g3 = adata[(adata.obs['bmi_lv'] > 20) & (adata.obs['bmi_lv'] < 25), :]
g4 = adata[adata.obs['bmi_lv'] < 20, :]

def categorize_bmi(adata):
    adata.obs['bmi_groups'] = pd.cut(adata.obs['bmi_lv'], bins=[0, 20, 25, 30, np.inf], labels=['bmi_<20', 'bmi_20-25', 'bmi_25-30', 'bmi_30+'])
    return adata

groups = [g1, g2, g3, g4]

highest_expr_genes = log(sc.pl.highest_expr_genes)
obesity_genes = [tup[-1] for tup in adata.uns['rank_genes_groups']['names']][:20]
anti_obesity_genes = [tup[0] for tup in adata.uns['rank_genes_groups']['names']][:20]

for group, name in zip(groups, ['bmi_30+', 'bmi_25-30', 'bmi_20-25', 'bmi_<20']):
    print(f"Number of cells in group {group.obs['bmi_lv'].unique()[0]}: {len(group)}")
    # plot highest expressed genes
    highest_expr_genes(group, n_top=20, save=f'_{DATASET}_{name}.png')

    plot_celltypes(group)

    sc.pl.umap(group, color='RNA.Class.Jun21_2024', save=f'_{DATASET}_{name}_RNA.Class.Jun21_2024.png')
    # sc.pl.umap(group, color='RNA.Subtype.Jun21_2024', save=f'_{DATASET}_RNA.Subtype.Jun21_2024.png')
    # sc.pl.umap(group, color='RNA.Subclass.Jun21_2024', save=f'_{DATASET}_RNA.Subclass.Jun21_2024.png')
