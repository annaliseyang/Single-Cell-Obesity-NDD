import scanpy as sc
import os
from utils import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import json

###################################

DATASET = "AD427"
PATH = "/home/anna_y/data/HumanBrainObesity/snRNA/"
FILENAME = "rna2.AD427_ADMR.QC.3343094.Jun24_2024.h5ad"

PROCESSED_DATA_PATH = "/home/anna_y/data/processed/"

###################################

def init_settings():
    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor="white")

@log
def load_data(dataset, path):
    adata = sc.read(path, cache=True)
    # adata.var_names_make_unique()

    # results_file = f"/home/anna_y/data/write/{dataset}.h5ad"  # the file that will store the analysis results
    # adata.write(results_file)

    return adata

init_settings()
# adata = load_data(DATASET, f'/home/anna_y/data/write/processed_{DATASET}_Jun17_2024.h5ad')
adata = load_data(DATASET, '/home/anna_y/data/write/AD427_ADMR_meta_Jul22_2024.h5ad')
print("Data loaded!", flush=True)
print(adata, flush=True)
