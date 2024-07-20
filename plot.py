from init import *

# plotting function
# name, func = 'TSNE', log(sc.pl.tsne)
# name, func = 'PCA', log(sc.pl.pca)
# name, func = 'UMAP', log(sc.pl.umap)
name, func = 'heatmap', log(sc.pl.heatmap)

variables = ['leiden_clust', 'bmi_lv', 'log1p_total_counts', 'log1p_n_genes_by_counts', 'region', 'batch']
highest_expr_genes = ['MALAT1', 'PCDH9', 'KCNIP4', 'CADM2', 'IL1RAPL1', 'DLG2', 'МТ-СО3', 'МТ-СО2', 'RBFOX1', 'SNHG14', 'CSMD1', 'CNTNAP2', 'МТ-АТР6', 'MAGI2', 'NRXN3', 'PTPRD', 'LSAMP', 'ADGRB3', 'NRG3', 'MT-CO1']
highest_expr_genes = [gene for gene in highest_expr_genes if gene in adata.var_names]
# obesity_genes = ['AC098829.1', 'ZNF638', 'SIK3', 'FAM118A', 'RPL37A', 'TAOK3', 'MACF1', 'PER3', 'INTU', 'SPRED2', 'GNG7', 'PPP6R3', 'TAF1D', 'HDAC4', 'RERE', 'MYO6', 'SEPTIN2', 'MARK3', 'ARGLU1', 'UBR2']
# obesity_genes = [gene for gene in obesity_genes if gene in adata.var_names]
cell_types = ['RNA.Subtype.Jun21_2024', 'RNA.Class.Jun21_2024', 'RNA.Subclass.Jun21_2024']

def plot(adata, key):
    if f'{name.lower()}_{DATASET}_{key}.png' in os.listdir('/home/anna_y/scRNA/figures/'):
        print(f"{name} for {key} already exists!", flush=True)
        return

    print(f"Generating {name} for {key}...", flush=True)
    try:
        func(adata, color=key, vmax='p99', save=f'_{DATASET}_{key}.png')
        print(f"{name} saved for {key}!", flush=True)
    except Exception as e:
        print(f"Error saving {name} for {key}: {e}", flush=True)

def plot_each(adata, keys):
    for key in keys:
        if f'{name.lower()}_{DATASET}_{key}.png' in os.listdir('/home/anna_y/scRNA/figures/'):
            print(f"{name} for {key} already exists!", flush=True)
            continue

        print(f"Generating {name} for {key}...", flush=True)
        try:
            func(adata, color=key, vmax='p99', save=f'_{DATASET}_{key}.png')
            print(f"{name} saved for {key}!", flush=True)
        except Exception as e:
            print(f"Error saving {name} for {key}: {e}", flush=True)

def plot_all(adata, keys, save=None):
    keys = [key for key in keys if key in adata.var_names]
    try:
        func(adata, color=keys, vmax='p99', save=f'_{DATASET}_highest_expr_genes.png')
        print(f"{name} saved for highest_expr_genes!", flush=True)
    except Exception as e:
        print(f"Error saving {name} for highest_expr_genes: {e}", flush=True)


def plot_celltype(adata, celltype: str, key: str):
    subset = adata[adata.obs['RNA.Class.Jun21_2024'] == celltype, :]
    try:
        func(subset, color=key, save=f'_{DATASET}_{celltype}_{key}.png')
        print(f"{name} saved for {celltype}!", flush=True)
    except Exception as e:
        print(f"Error saving {name} for {celltype}: {e}", flush=True)

def plot_celltypes(adata):
    for type in adata.obs['RNA.Class.Jun21_2024'].unique():
        plot_celltype(adata, type, 'RNA.Subclass.Jun21_2024')


if __name__ == "__main__":
    # plot_all(adata, highest_expr_genes)
    # plot_celltypes(adata)
    sc.pl.dotplot(adata, highest_expr_genes, groupby='RNA.Class.Jun21_2024', standard_scale='var', save=f'{DATASET}_heg_RNA.Class.Jun21_2024.png')
    sc.pl.matrixplot(adata, highest_expr_genes, groupby='RNA.Class.Jun21_2024', standard_scale='var', save=f'{DATASET}_heg_RNA.Class.Jun21_2024.png')
