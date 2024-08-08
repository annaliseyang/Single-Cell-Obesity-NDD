library(SeuratDisk)

# path = "/home/anna_y/data/write/"
# name = "AD427_ADMR_Aug6_2024"
# filename = "AD427_ADMR_Aug6_2024.h5ad"

argv <- commandArgs(T)
filepath = argv[1]

# Convert .h5ad file to .h5seurat
Convert(filepath, dest = "h5seurat", overwrite = TRUE)
print(paste0(filepath, " file converted to .h5seurat"))
