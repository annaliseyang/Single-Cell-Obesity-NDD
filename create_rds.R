library(Seurat)
library(data.table)
library(readr)
library(Matrix)

dir = commandArgs(T)[1] # directory that contains the csv files from running export_info.py
files = list.files(dir, pattern=".csv")

# Load the dense matrix (gene expression data)
counts_filename = list.files(dir, pattern="expression_matrix.csv")[1]
counts_path = file.path(dir, counts_filename)
print(paste("Loading counts matrix from: ", counts_path))

# Use a sparse matrix representation for the counts
counts <- read_csv_chunked(counts_path, chunk_size=10000, col_types = cols(.default = col_double()), callback=DataFrameCallback$new(function(chunk, pos) {
  chunk_sparse <- Matrix(as.matrix(chunk), sparse = TRUE)
  print(dim(chunk_sparse))
  return(chunk_sparse)
}))

# Transpose the counts matrix
counts <- t(counts)

seurat_object <- CreateSeuratObject(counts = counts)
print(seurat_object)

# Add metadata
metadata_filename <- list.files(dir, pattern="obs.csv")[1]
print(paste("Loading metadata from:", metadata_filename))
meta_data <- read.csv(file.path(dir, metadata_filename), header = TRUE)
rownames(meta_data) <- meta_data$Sample_barcode

seurat_object <- AddMetaData(object = seurat_object, metadata = meta_data)

# Save the seurat object to an rds file
out_filename = sub("_expression_matrix.csv", ".rds", counts_filename)
out_path = file.path(dir, out_filename)
saveRDS(seurat_object, file = out_path)
print(paste0(".rds object saved to ", out_path))
