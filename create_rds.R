library(Seurat)
library(data.table)
library(readr)
library(Matrix)

# Input directory that contains the csv files from running export_info.py
dir <- commandArgs(trailingOnly = TRUE)[1]

# Load the dense matrix (gene expression data)
counts_filename <- list.files(dir, pattern="expression_matrix.csv")[1]
counts_path <- file.path(dir, counts_filename)
print(paste("Loading counts matrix from:", counts_path))

# Initialize an empty list to store chunks as sparse matrices
sparse_matrices <- list()

# Process the file in chunks and convert each chunk to a sparse matrix
read_csv_chunked(counts_path, chunk_size = 100000,
                 col_types = cols(.default = col_double()),
                 callback = DataFrameCallback$new(function(chunk, pos) {
  sparse_chunk <- Matrix(as.matrix(chunk), sparse = TRUE)
  sparse_matrices[[length(sparse_matrices) + 1]] <<- sparse_chunk
  print(paste("Processed chunk at position:", pos, "with dimensions:", dim(sparse_chunk)[1], dim(sparse_chunk)[2])
)
}))

# Combine all sparse matrices into one large sparse matrix
print("Combining sparse matrices...")
counts_sparse <- do.call(rbind, sparse_matrices)

# Transpose the counts matrix
print("Transposing sparse matrix...")
counts_sparse <- t(counts_sparse)

# Create Seurat object
seurat_object <- CreateSeuratObject(counts = counts_sparse)
print(seurat_object)

# Add metadata
metadata_filename <- list.files(dir, pattern="obs.csv")[1]
print(paste("Loading metadata from:", metadata_filename))
meta_data <- fread(file.path(dir, metadata_filename))
rownames(meta_data) <- meta_data$Sample_barcode

seurat_object <- AddMetaData(object = seurat_object, metadata = meta_data)

# Save the Seurat object to an RDS file
out_filename <- sub("_expression_matrix.csv", ".rds", counts_filename)
out_path <- file.path(dir, out_filename)
saveRDS(seurat_object, file = out_path)
print(paste0(".rds object saved to ", out_path))
