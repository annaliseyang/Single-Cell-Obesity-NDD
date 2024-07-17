import scanpy as sc
import os
from utils import *

###################################

DATASET = "AD427"
PATH = "/home/anna_y/data/HumanBrainObesity/snRNA/"
FILENAME = "rna2.AD427_ADMR.QC.3343094.Jun24_2024.h5ad"

###################################

def init_settings():
    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor="white")

@log
def load_data(dataset, path):
    adata = sc.read_h5ad(path)
    adata.var_names_make_unique()

    # results_file = f"/home/anna_y/data/write/{dataset}.h5ad"  # the file that will store the analysis results
    # adata.write(results_file)

    return adata

@log
def preprocess(adata):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 10, :]

    return adata

if __name__ == "__main__":
    init_settings()

    adata = load_data(DATASET, PATH + FILENAME)
    results_file = f"/home/anna_y/data/write/{DATASET}_new.h5ad"  # the file that will store the analysis results
    adata.write(results_file)

    print(adata, flush=True)

    adata = preprocess(adata)
    print(adata, flush=True)

    # # annotate the group of mitochondrial genes as "mt"
    # adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # sc.pp.calculate_qc_metrics(
    #     adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    # )
    sc.pp.normalize_total(adata, target_sum=1e4)

    # sc.pl.clustermap(adata)
    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", save=f'_{DATASET}_pct_counts_mt.png')
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", save=f'_{DATASET}_n_genes_by_counts.png')

    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    sc.pl.highly_variable_genes(adata, save=f'_{DATASET}_highly_variable_genes.png')
    adata.raw = adata
