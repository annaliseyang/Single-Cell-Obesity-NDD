library(Seurat)

dir = commandArgs(T)[1] # directory that contains the csv files from running export_info.py
files = list.files(dir, pattern=".csv")

# Load the dense matrix (gene expression data)
counts_filename = list.files(dir, pattern="expression_matrix.csv")[1]
counts_path = paste0(dir, counts_filename)
counts <- read.csv(counts_path, row.names = 1)

# # Load the cell metadata
# obs_filename = list.files(dir, pattern="obs.csv")[1]
# obs_path = paste0(dir, obs_filename)
# obs <- read.csv(obs_path, row.names = 1)

# # Load the gene metadata
# var_filename = list.files(dir, pattern="var.csv")[1]
# var_path = paste0(dir, obs_filename)
# var <- read.csv(var_path, row.names = 1)

# Create a Seurat object
# seurat_object <- CreateSeuratObject(counts = as.matrix(counts), meta.data = obs)
counts <- t(counts)
seurat_object <- CreateSeuratObject(counts = as.matrix(counts))
print(seurat_object)

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
