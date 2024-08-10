library(Seurat)
library(data.table)
library(readr)

dir = commandArgs(T)[1] # directory that contains the csv files from running export_info.py
files = list.files(dir, pattern=".csv")

# Load the dense matrix (gene expression data)
counts_filename = list.files(dir, pattern="expression_matrix.csv")[1]
counts_path = paste0(dir, counts_filename)
print(paste("Loading counts matrix from: ", counts_path))
# counts <- read.csv(counts_path, row.names = 1)

process_chunk <- function(chunk, pos) {
  # Do some processing with the chunk
  print(dim(chunk))
  # Return the chunk or process it
  return(chunk)
}

counts <- read_csv_chunked(counts_path, chunk_size=10000, callback=DataFrameCallback$new(process_chunk))
summary(counts)
head(counts)

counts <- t(counts) # counts matrix column names must be cell names
seurat_object <- CreateSeuratObject(counts = as.matrix(counts))
print(seurat_object)

# Add metadata
metadata_filename <- list.files(dir, pattern="obs.csv")[1]
print(paste("Loading metadata from:", metadata_filename))
meta_data <- read.csv(file.path(dir, metadata_filename), header = TRUE)
rownames(meta_data) <- meta_data$Sample_barcode

seurat_object <- AddMetaData(object = seurat_object, metadata = meta_data)

# Save the seurat object to an rds file
out_filename = sub("_expression_matrix.csv", ".rds", counts_filename)
out_path = paste0(dir, out_filename)
saveRDS(seurat_object, file = out_path)
print(paste0(".rds object saved to ", out_path))

# cell_names = rownames(seurat_object)
# rownames(obs) <- cell_names


# Add gene metadata to the Seurat object
# seurat_object@meta.data <- obs
# seurat_object@assays$RNA@meta.features <- var
