import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import os
import sys

def cell_fraction_discrete(adata, groupby='bmi_groups', celltype_col='Class', total=None, save=None):
    cell_types = adata.obs[celltype_col]
    groups = adata.obs[groupby]

    # Calculate fractions
    cell_fractions = adata.obs.groupby([groups, cell_types]).size()
    cell_fractions = cell_fractions.unstack(fill_value=0)
    total = total if total else cell_fractions.sum(axis=1)
    cell_fractions = cell_fractions.div(total, axis=0)
    print(cell_fractions.head())

    # # Calculate standard deviation for error bars
    # cell_fractions_std = adata.obs.groupby([groups, cell_types]).size().unstack(fill_value=0).div(
    #     adata.obs.groupby(groups).size(), axis=0
    # ).groupby(cell_types).std()
    # print(cell_fractions_std.head())

    # # Convert to long format for seaborn
    # cell_fractions_long = cell_fractions.reset_index().melt(id_vars=groups)
    # cell_fractions_std_long = cell_fractions_std.reset_index().melt(id_vars=groups)

    # Plot
    celltypes_unique = adata.obs[celltype_col].unique()
    width, height = 3, 4
    nrow, ncol = 1, len(celltypes_unique)
    fig, axs = plt.subplots(nrows=nrow, ncols=ncol, figsize=(ncol*width, nrow*height), sharex=True, sharey=True)
    for i, col in enumerate(celltypes_unique):
        subset = adata[adata.obs[celltype_col] == col]
        sns.barplot(data=cell_fractions, x=groupby, y=col, hue=col, ax=axs[i])
        axs[i].set_ylabel(f'{col} Fraction')
        axs[i].set_title(f'{col}')
        axs[i].tick_params(axis='x', rotation=45)
        axs[i].legend_.remove()
    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=300)
        print(f"Cell fraction saved to {save}")


def compute_cell_fractions(adata, total=None, celltype_col='Class'):
    total_cells = total if total else adata.obs.shape[0]
    cell_fractions = adata.obs.groupby(celltype_col, observed=False).size() / total_cells
    return cell_fractions

# def cell_fraction_discrete(adata, groupby='bmi_groups_0', celltype_col='Class', save=None):
#     celltypes = adata.obs[celltype_col].unique()
#     groups = adata.obs[groupby].unique()
#     for group in groups:
#         subset = adata[adata.obs[groupby] == group, :]
#         cell_fractions = compute_cell_fractions(subset, celltype_col=celltype_col)

def cell_fraction_continuous(adata, x='bmi_lv', celltype_col='Class', total=None, save=None):
    celltypes = adata.obs[celltype_col].unique()
    x_values = adata.obs[x].unique()
    projids = adata.obs['projid'].unique()

    df = pd.DataFrame(index=projids, columns=[x] + list(celltypes))

    for projid in projids:
        subset = adata[adata.obs['projid'] == projid, :]
        x_val = subset.obs[x].iloc[0]
        total_cells = subset.shape[0] # total cells for this individual
        print(f"Processing projid {projid} with {x}: {x_val}, total cells: {total_cells=}")

        cell_fractions = compute_cell_fractions(subset, celltype_col=celltype_col, total=total_cells)
        # print(f"Cell fractions for {projid}: {cell_fractions}")
        df.at[projid, x] = x_val
        df.loc[projid, celltypes] = cell_fractions
        # break
    print(df.head())

    # scatter plot
    width, height = 3, 3
    ncol = 7
    nrow = int(np.ceil(len(celltypes) / ncol))

    fig, axs = plt.subplots(nrows=nrow, ncols=ncol, figsize=(ncol*width, nrow*height), sharex=True, sharey=True)
    # celltype='Oli'
    for i, celltype in enumerate(celltypes):
        row, col = divmod(i, ncol)
        ax_index = (row, col) if nrow > 1 else i
        sns.scatterplot(data=df, x=x, y=celltype, ax=axs[ax_index])

        x_values = df[x].astype(float)
        y_values = df[celltype].astype(float)
        valid_indices = ~np.isnan(x_values) & ~np.isnan(y_values)
        x_values = x_values[valid_indices]
        y_values = y_values[valid_indices]

        sns.regplot(x=x_values, y=y_values, ax=axs[ax_index], scatter=False)

        if len(x_values.dropna()) > 1 and len(y_values.dropna()) > 1:
            corr, p_value = pearsonr(x_values, y_values)
            print(f"Correlation between {x} and {celltype}: corr={corr}, p-value={p_value}")
            p_label = ('*' if p_value < 0.05 else '') + f"p={round(p_value, 3)}"
            # axs[ax_index].annotate(p_label, xy=(0.05, 0.95), xycoords='axes fraction',
            #                        fontsize=10, color='black', verticalalignment='top')

        axs[ax_index].set_xlabel(x)
        axs[ax_index].set_ylabel('Fraction')
        axs[ax_index].set_title(f'{celltype} ({p_label})')
    # fig.set_title(f"Cell fractions vs. {x}")
    plt.tight_layout()

    if save:
        if not os.path.exists(os.path.dirname(save)):
            os.makedirs(os.path.dirname(save))
        plt.savefig(save, dpi=300)
        print(f"Cell fraction scatterplot saved to {save}.")


if __name__ == "__main__":
    # in_path = "/home/anna_y/data/test/tiny_AD427MR/tiny_AD427MR.h5ad"
    # in_path = "/home/anna_y/data/write/all/AD427MR_50k/AD427MR_50k.h5ad"
    in_path = "/home/anna_y/data/write/all/AD427MR/AD427MR.h5ad"

    # in_path = sys.argv[1]  # e.g. /home/anna_y/data/write/all/AD427MR/AD427MR.h5ad
    adata = sc.read_h5ad(in_path)
    print(adata)

    x = 'bmi_lv'
    # x = 'gpath'
    adata = adata[adata.obs[x].notna(), :].copy() # filter cells

    for celltype_col in ['Class', 'Subclass', 'Subtype']:
        out_path = in_path.replace('write', 'results/figures/cell_fraction').replace('.h5ad', f'_cell_fraction_{x}_{celltype_col}_cont.png')
        print("out_path:", out_path)
        cell_fraction_continuous(adata, x=x, celltype_col=celltype_col, save=out_path)
