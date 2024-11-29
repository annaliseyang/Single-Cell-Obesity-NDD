indir = commandArgs(T)[1] # eg. /home/anna_y/data/write/Class/Ast/
filename <- list.files(indir, pattern="_bmi.rds")[1]
region = "PFC"
out_filename <- sub("_bmi.rds", paste0('.', region, '.rds'), filename)
out_dir <- file.path(indir, out_filename)

print(paste("Reading seurat object from:", filename))
seurat_obj <- readRDS(file.path(indir, filename))
metadata <- seurat_obj@meta.data

# filter cells by batch
print(unique(metadata$region))
filtered_cells <- rownames(metadata[metadata$region=="PFC", ])
# print("metadata$batch:", metadata$batch)
print(paste0("Number of cells after filtering by region=='PFC': ", length(filtered_cells)))

# create Seurat object with filtered cells
seurat_obj_filtered <- subset(seurat_obj, cells = filtered_cells)

# save
print(paste0("Saved filtered Seurat object to: ", out_dir))
saveRDS(seurat_obj_filtered, out_dir)
print("Finished.")
