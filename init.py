import scanpy as sc
import os
from utils import *
from preprocessing import *
import numpy as np

init_settings()
adata = load_data(DATASET, PATH + FILENAME)
