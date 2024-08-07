library(SeuratDisk)

argv <- commandArgs(T)
filepath = argv[1]
out_path = gsub(".h5seurat", ".rds", filepath)

# Save the Seurat object as an RDS file
seurat_object <- LoadH5Seurat(paste0(filepath, ".h5seurat"))
saveRDS(seurat_object, file = out_path)
print(paste0("Seurat object saved as ", name, ".rds"))
