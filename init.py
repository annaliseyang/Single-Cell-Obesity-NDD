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

DATASET = "All"
PATH = "/home/anna_y/data/HumanBrainObesity/snRNA/"
FILENAME = "rna2.AD427_ADMR.QC.3343094.Jun24_2024.h5ad"

PROCESSED_DATA_PATH = "/home/anna_y/data/processed/"
in_path = sys.argv[1]

###################################

def init_settings():
    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor="white")

init_settings()
adata = sc.read_h5ad(in_path)
print("Data loaded!", flush=True)
print(adata, flush=True)
