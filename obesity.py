from plot import *
from deg import get_top_n_genes

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
AD_groups = ['earlyAD', 'lateAD', 'nonAD']
AD_subsets = {state: adata[adata.obs['ADdiag3types'].isin([state]) | adata.obs['Pathology'].isin([state])] for state in AD_groups}

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
obesity_genes = get_top_n_genes('bmi_30+', 'rank_genes_groups_bmi.csv', n_genes=20)
anti_obesity_genes = get_top_n_genes('bmi_<20', 'rank_genes_groups_bmi.csv', n_genes=20)
early_AD_genes = get_top_n_genes('earlyAD', 'rank_genes_groups_AD.csv', n_genes=20)
late_AD_genes = get_top_n_genes('lateAD', 'rank_genes_groups_AD.csv', n_genes=20)

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
def umap_by_groups(subsets, color='RNA.Class.Jun21_2024', groupby:str=None, save=None):
    side_length = 4
    fig, axs = plt.subplots(1, len(subsets), figsize=(len(subsets) * side_length, side_length), sharex=True, sharey=True)
    vmax = np.percentile(np.concatenate([subset[:, color].X.toarray().flatten() for subset in subsets.values()]), 99)
    for i, group in enumerate(subsets.keys()):
        subset = subsets[group]
        sc.pl.umap(subset, color=color, vmax=vmax, ax=axs[i], show=False)
        axs[i].set_title(f"{color}, {group}")
    plt.tight_layout()
    if save:
        plt.savefig(save)
    else:
        plt.savefig(f"figures/{groupby}/umap_{color}_by_{groupby}.png")

@log
def umap_by_bmi(bmi_subsets, color='RNA.Class.Jun21_2024', save=None):
    fig, axs = plt.subplots(1, 4, figsize=(20, 5), sharex=True, sharey=True)
    vmax = np.percentile(np.concatenate([subset[:, color].X.toarray().flatten() for subset in bmi_subsets.values()]), 99)

    for i, group in enumerate(bmi_groups):
        subset = bmi_subsets[group]
        sc.pl.umap(subset, color=color, vmax=vmax, ax=axs[i], show=False)
        axs[i].set_title(f"{color}, {group}")

    plt.tight_layout()
    if save:
        plt.savefig(save)
    else:
        plt.savefig(f"figures/bmi/umap_{color}_p99_by_bmi.png")


def umap_top_n_genes_by_groups(subsets, colors, n_top=5, groupby:str=None, save=None):
    side_length = 4
    fig, axs = plt.subplots(n_top, len(subsets), figsize=(len(subsets) * side_length, side_length*n_top), sharex=True, sharey=True)
    for row, color in enumerate(colors[:n_top]):
        vmax = np.percentile(np.concatenate([subset[:, color].X.toarray().flatten() for subset in subsets.values()]), 99)
        for col, group in enumerate(subsets.keys()):
            subset = bmi_subsets[group]
            try:
                sc.pl.umap(subset, color=color, vmax=vmax, ax=axs[row, col], show=False)
                # sc.pl.dotplot(subset, obesity_genes[:n_top], groupby='bmi_groups', standard_scale='var', ax=axs[row, col])
                axs[row, col].set_title(f"{color}, {group}")
            except Exception as e:
                print(f"Error plotting {color} for {group}: {e}")

    plt.tight_layout()
    if save:
        plt.savefig(save)
    else:
        plt.savefig(f"figures/{groupby}/umap_top_{n_top}_genes_by_{groupby}.png")

@log
def umap_top_n_genes_by_bmi(bmi_subsets, colors, n_top=5, save=None):
    fig, axs = plt.subplots(n_top, len(bmi_groups), figsize=(20, 5*n_top), sharex=True, sharey=True)
    for row, color in enumerate(colors[:n_top]):
        vmax = np.percentile(np.concatenate([subset[:, color].X.toarray().flatten() for subset in bmi_subsets.values()]), 99)
        for col, group in enumerate(bmi_groups):
            subset = bmi_subsets[group]
            try:
                sc.pl.umap(subset, color=color, vmax=vmax, ax=axs[row, col], show=False)
                # sc.pl.dotplot(subset, obesity_genes[:n_top], groupby='bmi_groups', standard_scale='var', ax=axs[row, col])
                axs[row, col].set_title(f"{color}, {group}")
            except Exception as e:
                print(f"Error plotting {color} for {group}: {e}")

    plt.tight_layout()
    if save:
        plt.savefig(save)
    else:
        plt.savefig(f"figures/bmi/umap_top_{n_top}_genes_by_bmi.png")


@log
def heatmap_by_bmi(adata, n_genes=5):
    groups = adata.obs['RNA.Class.Jun21_2024'].unique()
    fig, axs = plt.subplots(ncols=len(groups), sharex=True, sharey=True)
    for i, group in enumerate(groups):
        subset = adata[adata.obs['RNA.Class.Jun21_2024'] == group, :]
        try:
            sc.pl.rank_genes_groups_heatmap(subset, n_genes=n_genes, groupby='RNA.Subclass.Jun21_2024', standard_scale='var', ax=axs[i], save=f'_obesity_DE_{group}_RNA.Subclass.Jun21_2024.png')
            axs[1, i].set_title(f"{group}")
        except Exception as e:
            print(f"An error occured: {e}")
            sc.pl.rank_genes_groups_heatmap(subset, n_genes=n_genes, groupby='RNA.Subclass.Jun21_2024', standard_scale='var', save=f'_obesity_DE_{group}_RNA.Subclass.Jun21_2024.png')


def deg_by_AD_state(subsets, groupby='AD_states', reference='nonAD'):
    # Perform DEG analysis using the rank_genes_groups function
    for group, subset in subsets.items():
        sc.tl.rank_genes_groups(subset, groupby=groupby, reference=reference)
    # Access the DEG results for each AD group
    for group, subset in subsets.items():
        print(f"DEG results for {group}:")
        print(subset.uns['rank_genes_groups'])


if __name__ == "__main__":
    # adata_bmi = adata[adata.obs['bmi_lv'].notna(), :]
    # sc.pl.dotplot(adata, obesity_genes, groupby='bmi_groups', standard_scale='var', save=f'{DATASET}_obesity_genes_bmi_groups.png')
    # sc.pl.dotplot(adata, anti_obesity_genes, groupby='bmi_groups', standard_scale='var', save=f'{DATASET}_anti_obesity_genes_bmi_groups.png')
    # sc.pl.dotplot(adata, obesity_genes, groupby='RNA.Class.Jun21_2024', standard_scale='var', save=f'{DATASET}_obesity_genes_RNA.Class.Jun21_2024.png')

    # sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, groupby='RNA.Class.Jun21_2024', standard_scale='var', save='_obesity_DE_RNA.Class.Jun21_2024.png')
    # sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, groupby='RNA.Subclass.Jun21_2024', standard_scale='var', save='_obesity_DE_RNA.Subclass.Jun21_2024.png')
    # sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, groupby='RNA.Subtype.Jun21_2024', standard_scale='var', save='_obesity_DE_RNA.Subtype.Jun21_2024.png')
    # sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, groupby='bmi_groups', standard_scale='var', save='_obesity_DE_bmi_groups.png')

    # for gene in obesity_genes:
    #     umap_by_bmi(bmi_subsets, color=gene)
    # for gene in late_AD_genes:
    #     umap_by_groups(AD_subsets, color=gene, save=f"figures/AD/umap_{gene}.png")

    # umap_top_n_genes_by_bmi(bmi_subsets, colors=obesity_genes, n_top=10)
    # umap_top_n_genes_by_bmi(bmi_subsets, colors=anti_obesity_genes, n_top=5, save=f"figures/bmi/umap_antiobesity_genes_by_bmi.png")
    # umap_top_n_genes_by_groups(AD_subsets, colors=early_AD_genes, groupby='AD', n_top=5)
    umap_top_n_genes_by_groups(AD_subsets, colors=late_AD_genes, groupby='AD', n_top=5)
    # heatmap_by_bmi(adata, n_genes=5)

    pass
