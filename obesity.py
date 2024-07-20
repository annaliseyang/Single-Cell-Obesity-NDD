from plot import *

# if __name__ == "__main__":
# adata = load_data(f'/home/anna_y/data/write/processed_{DATASET}_Jun17_2024.h5ad')
# g1 = adata[adata.obs['bmi_lv'] > 30, :]
# g2 = adata[(adata.obs['bmi_lv'] > 25) & (adata.obs['bmi_lv'] < 30), :]
# g3 = adata[(adata.obs['bmi_lv'] > 20) & (adata.obs['bmi_lv'] < 25), :]
# g4 = adata[adata.obs['bmi_lv'] < 20, :]

def categorize_bmi(adata, num_bins=4):
    adata.obs['bmi_groups'] = pd.cut(adata.obs['bmi_lv'], bins=[0, 20, 25, 30, np.inf], labels=['bmi_<20', 'bmi_20-25', 'bmi_25-30', 'bmi_30+'])
    return adata

# groups = [g1, g2, g3, g4]
bmi_groups = ['bmi_<20', 'bmi_20-25', 'bmi_25-30', 'bmi_30+']
bmi_subsets = {group: adata[adata.obs['bmi_groups'] == group] for group in bmi_groups}


def print_cell_counts_and_percentages(adata):
    for group in bmi_groups:
        # print(f"Number of people with {group}", len(adata[adata.obs['bmi_groups'] == group]))
        num_cells = len(adata[adata.obs['bmi_groups'] == group])
        print(f"Number of cells {group}:", num_cells)
        print(f"Percentage of cells {group}: {num_cells / len(adata) * 100:.2f}%\n")

    for state in ['earlyAD', 'lateAD', 'nonAD']:
        num_cells = len(adata[adata.obs['ADdiag3types'].isin([state]) | adata.obs['Pathology'].isin([state])])
        print(f"Number of cells {state}:", num_cells)
        print(f"Percentage of cells {state}: {num_cells / len(adata) * 100:.2f}%\n")

highest_expr_genes = log(sc.pl.highest_expr_genes)
obesity_genes = [tup[-1] for tup in adata.uns['rank_genes_groups']['names']][:20]
anti_obesity_genes = [tup[0] for tup in adata.uns['rank_genes_groups']['names']][:20]

def plot_celltypes_by_bmi(adata):
    for group in bmi_groups:
        subset = adata[adata.obs['bmi_groups'] == group]
        sc.pl.umap(subset, color='RNA.Class.Jun21_2024', save=f'_{group}_RNA.Class.Jun21_2024.png')
        # subset = adata[adata.obs['bmi_groups'] == group]
        # sc.pl.umap(subset, obesity_genes, save=f'_{group}_obesity_genes.png')

def highest_expr_genes_by_bmi(adata, n_top=20, save=None):
    for name, group in bmi_subsets.items():
        print(f"Number of cells in group {group.obs['bmi_lv'].unique()[0]}: {len(group)}")
        # plot highest expressed genes
        highest_expr_genes(group, n_top=20, save=f'_{DATASET}_{name}.png')

        plot_celltypes(group)

        sc.pl.umap(group, color='RNA.Class.Jun21_2024', save=f'_{DATASET}_{name}_RNA.Class.Jun21_2024.png')
        # sc.pl.umap(group, color='RNA.Subtype.Jun21_2024', save=f'_{DATASET}_RNA.Subtype.Jun21_2024.png')
        # sc.pl.umap(group, color='RNA.Subclass.Jun21_2024', save=f'_{DATASET}_RNA.Subclass.Jun21_2024.png')

@log
def umap_by_bmi(bmi_subsets, color='RNA.Class.Jun21_2024'):
    fig, axs = plt.subplots(1, 4, figsize=(20, 5), sharex=True, sharey=True)
    # vmin = min(adata.obs[color])
    # vmax = max(adata.obs[color]) * 0.99

    for i, group in enumerate(bmi_groups):
        subset = bmi_subsets[group]
        sc.pl.umap(subset, color=color, vmax='p99', ax=axs[i], show=False)
        axs[i].set_title(f"{color}, {group}")

    plt.tight_layout()
    plt.savefig(f"figures/bmi/umap_{color}_by_bmi.png")
    plt.show()

def umap_top_n_genes_by_bmi(bmi_subsets, colors, n_top=5):
    fig, axs = plt.subplots(n_top, 4, figsize=(20, 5*n_top), sharex=True, sharey=True)
    for row, color in enumerate(colors[:n_top]):
        # vmin = min(adata.obs[color])
        # vmax = max(adata.obs[color]) * 0.99

        for col, group in enumerate(bmi_groups):
            subset = bmi_subsets[group]
            sc.pl.umap(subset, color=color, vmax='p99', ax=axs[row, col], show=False)
            # sc.pl.dotplot(subset, obesity_genes[:n_top], groupby='bmi_groups', standard_scale='var', ax=axs[row, col])
            axs[row, col].set_title(f"{color}, {group}")

        plt.tight_layout()
        plt.savefig(f"figures/bmi/umap_top_{n_top}_genes_by_bmi.png")

if __name__ == "__main__":
    # adata_bmi = adata[adata.obs['bmi_lv'].notna(), :]
    # sc.pl.dotplot(adata, obesity_genes, groupby='bmi_groups', standard_scale='var', save=f'{DATASET}_obesity_genes_bmi_groups.png')
    # sc.pl.dotplot(adata, anti_obesity_genes, groupby='bmi_groups', standard_scale='var', save=f'{DATASET}_anti_obesity_genes_bmi_groups.png')
    # sc.pl.dotplot(adata, obesity_genes, groupby='RNA.Class.Jun21_2024', standard_scale='var', save=f'{DATASET}_obesity_genes_RNA.Class.Jun21_2024.png')

    # sc.pl.rank_genes_groups_heatmap(adata_bmi, n_genes=5, groupby='RNA.Class.Jun21_2024', standard_scale='var', save='_obesity_DE_RNA.Class.Jun21_2024.png')
    # sc.pl.rank_genes_groups_heatmap(adata_bmi, n_genes=5, groupby='bmi_groups', standard_scale='var', save='_obesity_DE_bmi_groups.png')

    for gene in obesity_genes:
        umap_by_bmi(bmi_subsets, color=gene)
    # gene = obesity_genes[sys.argv[1]]
    # print(f"Plotting UMAP for gene {gene}", flush=True)
    # umap_by_bmi(bmi_subsets, color=gene)
    # umap_top_n_genes_by_bmi(bmi_subsets, colors=obesity_genes, n_top=5)

    # umap_by_bmi(bmi_subsets, color='RNA.Class.Jun21_2024')
