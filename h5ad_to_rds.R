# install.packages("reticulate")
library(reticulate)
library(Seurat)
library(anndata)
library(SeuratDisk)

in_path = "/home/anna_y/data/test/tiny_AD427_ADMR_Aug6_2024.h5ad"
out_path = sub(".h5ad", ".rds", in_path) # output file will be saved in the same directory as the input
# filename <- basename(filepath)

# Load the Python environment
use_python("~/.conda/envs/sc2024/bin/python")

# Load h5ad data
anndata <- import("anndata", convert = FALSE)
adata <- anndata$read_h5ad(in_path)
print(adata)

# Convert the AnnData object to a Seurat object
# seurat_object <- Convert(source = in_path, dest = "seurat")
# seurat_object <- CreateSeuratObject(counts = t(adata$X), project = "tiny_AD427_ADMR_Aug6_2024", min.cells = 3, min.features = 200)
# # seurat_object <- CreateSeuratObject(counts = as.matrix(adata$X), project = "tiny_AD427_ADMR_Aug6_2024")
# seurat_object <- CreateSeuratObject(counts = adata$X[1], project = "tiny_AD427_ADMR_Aug6_2024")


seurat_object <- as.Seurat(adata)
print(seurat_object)

# Save the Seurat object as an RDS file
saveRDS(seurat_object, file = out_path)
