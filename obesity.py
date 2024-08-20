from deg import get_top_n_genes, rank_genes_groups
from utils import log
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

def categorize_bmi(adata, num_bins=4):
    adata.obs['bmi_groups'] = pd.cut(adata.obs['bmi_lv'], bins=[0, 20, 25, 30, np.inf], labels=['bmi_<20', 'bmi_20-25', 'bmi_25-30', 'bmi_30+'])
    return adata


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


def plot_celltypes_by_bmi(adata):
    for group in bmi_groups:
        subset = adata[adata.obs['bmi_groups'] == group]
        sc.pl.umap(subset, color='RNA.Class.Jun21_2024', save=f'_{group}_RNA.Class.Jun21_2024.png')
        # subset = adata[adata.obs['bmi_groups'] == group]
        # sc.pl.umap(subset, obesity_genes, save=f'_{group}_obesity_genes.png')

def highest_expr_genes_by_bmi(adata, n_top=20, save=None):
    from plot import plot_celltypes
    for name, group in bmi_subsets.items():
        print(f"Number of cells in group {group.obs['bmi_lv'].unique()[0]}: {len(group)}")
        # plot highest expressed genes
        highest_expr_genes(group, n_top=20, save=f'_{name}.png')

        plot_celltypes(group)

        sc.pl.umap(group, color='RNA.Class.Jun21_2024', save=f'_{name}_RNA.Class.Jun21_2024.png')

@log
def umap_by_groups(subsets, color='RNA.Class.Jun21_2024', groupby:str=None, save=None):
    side_length = 4
    fig, axs = plt.subplots(1, len(subsets), figsize=(len(subsets) * side_length, side_length), sharex=True, sharey=True)
    vmax = np.percentile(np.concatenate([subset[:, color].X.toarray().flatten() for subset in subsets.values()]), 99)
    vmin = 0
    for i, group in enumerate(subsets.keys()):
        subset = subsets[group]
        sc.pl.umap(subset, color=color, vmax=vmax, vmin=vmin, ax=axs[i], show=False)
        axs[i].set_title(f"{color}, {group}")
    plt.tight_layout()
    if save:
        plt.savefig(save)
    else:
        plt.savefig(f"/home/anna_y/data/results/figures/{groupby}/umap_{color}_by_{groupby}.png")

@log
def umap_by_bmi(bmi_subsets, color='RNA.Class.Jun21_2024', save=None):
    fig, axs = plt.subplots(1, 4, figsize=(20, 5), sharex=True, sharey=True)
    vmax = np.percentile(np.concatenate([subset[:, color].X.toarray().flatten() for subset in bmi_subsets.values()]), 99)
    vmin = 0

    for i, group in enumerate(bmi_groups):
        subset = bmi_subsets[group]
        sc.pl.umap(subset, color=color, vmax=vmax, vmin=vmin, ax=axs[i], show=False)
        axs[i].set_title(f"{color}, {group}")

    plt.tight_layout()
    if save:
        plt.savefig(save)
    else:
        plt.savefig(f"/home/anna_y/data/results/figures/bmi/umap_{color}_p99_by_bmi.png")

@log
def umap_top_n_genes_by_groups(subsets, colors: list, n_top=5, groupby:str=None, save=None):
    width, height = 4, 3.5
    fig, axs = plt.subplots(n_top, len(subsets), figsize=(len(subsets) * width, n_top * height), sharex=True, sharey=True)
    for row, color in enumerate(colors[:n_top]):
        values = np.concatenate([subset[:, color].X.toarray().flatten() for subset in subsets.values()])
        vmax = np.percentile(values, 99)
        if vmax == 0:
            # If no expression at the 99th percentile, set vmax to max expression
            vmax = max(values)
        vmin = 0
        print(f"Plotting {color}... vmax: {vmax}, vmin: {vmin}")
        for col, group in enumerate(subsets.keys()):
            subset = subsets[group]
            try:
                sc.pl.umap(subset, color=color, vmax=vmax, vmin=vmin, ax=axs[row, col], show=False)
                axs[row, col].set_title(f"{color}, {group}")
            except Exception as e:
                print(f"Error plotting {color} for {group}: {e}")

    plt.tight_layout()
    if save:
        plt.savefig(save)
    else:
        plt.savefig(f"/home/anna_y/data/results/figures/{groupby}/umap_top_{n_top}_genes_by_{groupby}.png")

@log
def umap_top_n_genes_by_bmi(bmi_subsets, colors, n_top=5, save=None):
    umap_top_n_genes_by_groups(bmi_subsets, colors, n_top=n_top, groupby='bmi', save=save)

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


if __name__ == "__main__":
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
    #     umap_by_groups(AD_subsets, color=gene, save=f"/home/anna_y/data/results/figures/AD/umap_{gene}.png")

    # umap_top_n_genes_by_bmi(bmi_subsets, colors=obesity_genes, n_top=10)
    # umap_top_n_genes_by_bmi(bmi_subsets, colors=anti_obesity_genes, n_top=5, save=f"/home/anna_y/data/results/figures/bmi/umap_antiobesity_genes_by_bmi.png")
    # umap_top_n_genes_by_groups(AD_subsets, colors=early_AD_genes, groupby='AD', n_top=5, save=f"/home/anna_y/data/results/figures/AD/umap_early_AD_genes_by_AD_states.png")
    # umap_top_n_genes_by_groups(AD_subsets, colors=late_AD_genes, groupby='AD', n_top=5, save=f"/home/anna_y/data/results/figures/AD/umap_late_AD_genes_by_AD_states.png")

    in_path = sys.argv[1]
    adata = sc.read_h5ad(in_path)

    # groups = [g1, g2, g3, g4]
    bmi_groups = ['bmi_<20', 'bmi_20-25', 'bmi_25-30', 'bmi_30+']
    bmi_subsets = {group: adata[adata.obs['bmi_groups'] == group] for group in bmi_groups}
    AD_groups = ['earlyAD', 'lateAD', 'nonAD']
    AD_subsets = {state: adata[adata.obs['ADdiag3types'].isin([state]) | adata.obs['Pathology'].isin([state])] for state in AD_groups}
    celltypes = ['Exc', 'Inh', 'Oli', 'Ast', 'Mic_Immune', 'OPC', 'Vasc_Epithelia']
    celltype_subsets = {celltype: adata[adata.obs['RNA.Class.Jun21_2024'] == celltype, :] for celltype in celltypes}

    highest_expr_genes = log(sc.pl.highest_expr_genes)
    obesity_genes = get_top_n_genes('bmi_30+', 'rank_genes_groups_bmi.csv', n_genes=20)
    anti_obesity_genes = get_top_n_genes('bmi_<20', 'rank_genes_groups_bmi.csv', n_genes=20)
    early_AD_genes = get_top_n_genes('earlyAD', 'rank_genes_groups_AD.csv', n_genes=20)
    late_AD_genes = get_top_n_genes('lateAD', 'rank_genes_groups_AD.csv', n_genes=20)

    for type, subset in celltype_subsets.items():
        # make subsets for each cell type
        AD_subsets = {group: subset[subset.obs['AD_states'] == group, :] for group in AD_groups}
        bmi_subsets = {group: subset[subset.obs['bmi_groups'] == group, :] for group in bmi_groups}

        # plot umaps for the top 20 cell type specific AD and bmi markers
        # for group in ['earlyAD', 'lateAD']:
        #     markers = get_top_n_genes('earlyAD', f'rank_genes_groups_{type}_by_AD_states.csv', n_genes=20)
        #     print(f"Top 20 early AD genes for {type}: {markers}", flush=True)
        #     umap_top_n_genes_by_groups(AD_subsets, colors=markers, n_top=5, save=f'/home/anna_y/data/results/figures/celltypes/{type}_early_AD_genes_celltype_specific.png')

        groupby = 'bmi_groups'
        print(f"Computing DEG results for {type} by {groupby}...", flush=True)
        # sc.tl.rank_genes_groups(subset, groupby=groupby, method='wilcoxon', key_added=f'rank_genes_groups_{type}_by_{groupby}')
        rank_genes_groups(subset, groupby=groupby, type=type)

        for group in bmi_groups:
            # try:
            #     markers = get_top_n_genes(group, f'rank_genes_groups_{type}_by_bmi_groups.csv', n_genes=20)
            #     print(f"Top 20 bmi genes for {type}: {markers}", flush=True)
            #     umap_top_n_genes_by_groups(bmi_subsets, colors=markers, n_top=5, save=f'/home/anna_y/data/results/figures/celltypes/{type}_bmi_genes_celltype_specific.png')
            # except Exception as e:

            markers = get_top_n_genes(group, f'rank_genes_groups_{type}_by_{groupby}.csv', n_genes=20)
            umap_top_n_genes_by_groups(bmi_subsets, colors=markers, n_top=5, save=f'/home/anna_y/data/results/figures/celltypes/{type}_{group}_genes_celltype_specific.png')




        # # umap_by_groups(AD_subsets, 'SLC38A2', groupby='celltype', save=f'/home/anna_y/data/results/figures/celltypes/{type}_early_AD_SLC38A2.png')
        # for gene in early_AD_genes[:5]:
        #     umap_by_groups(AD_subsets, color=gene, groupby='celltype', save=f'/home/anna_y/data/results/figures/celltypes/{type}_early_AD_{gene}.png')
        # for gene in late_AD_genes[:5]:
        #     umap_by_groups(AD_subsets, color=gene, groupby='celltype', save=f'/home/anna_y/data/results/figures/celltypes/{type}_late_AD_{gene}.png')
        # for gene in obesity_genes[:5]:
        #     umap_by_groups(AD_subsets, color=gene, groupby='celltype', save=f'/home/anna_y/data/results/figures/celltypes/{type}_obesity_{gene}.png')
        # for gene in anti_obesity_genes[:5]:
        #     umap_by_groups(AD_subsets, color=gene, groupby='celltype', save=f'/home/anna_y/data/results/figures/celltypes/{type}_anti_obesity_{gene}.png')


    # heatmap_by_bmi(adata, n_genes=5)

    pass
