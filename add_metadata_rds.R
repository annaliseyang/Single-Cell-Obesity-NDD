library(Seurat)
print("Hello World!!")
print("How are you Anna!")
# Define file paths
# in_path <- "/home/anna_y/data/test/tiny_AD427_ADMR/"
in_path <- commandArgs(T)[1] # eg. /home/anna_y/data/write/Class/Ast/
rds_filename <- list.files(in_path, pattern=".rds")[1]
metadata_filename <- list.files(in_path, pattern="obs.csv")[1]

print(rds_filename)
print(metadata_filename)

# Load the Seurat object and metadata
seurat_object <- readRDS(file.path(in_path, rds_filename))
meta_data <- read.csv(file.path(in_path, metadata_filename), header = TRUE)

# Ensure the row names of metadata match the cell names in the Seurat object
rownames(meta_data) <- meta_data$Sample_barcode

# Check if all metadata row names are in the Seurat object cell names
if (!all(rownames(meta_data) %in% colnames(seurat_object))) {
    mismatched_cells <- setdiff(rownames(meta_data), colnames(seurat_object))
    stop("Mismatch between cell names in metadata and Seurat object. Please check the following cells: ", paste(mismatched_cells, collapse = ", "))
}

# Add metadata to the Seurat object
seurat_object <- AddMetaData(object = seurat_object, metadata = meta_data)

# View the updated Seurat object
cat("Number of cells:", ncol(seurat_object), "\n")
cat("Number of genes:", nrow(seurat_object), "\n")
head(seurat_object@meta.data)

# out_filename = sub(, ".rds", counts_filename)
out_filename = sub(".rds", "_meta.rds", rds_filename)
out_path = file.path(in_path, out_filename)
saveRDS(seurat_object, file = out_path)
print(paste("New .rds file with metadata saved to", out_path))
