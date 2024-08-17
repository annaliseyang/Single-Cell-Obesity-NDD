library(Seurat)

in_path = commandArgs(T)[1] # path to .rds file
# in_path = "~/data/write/Subtype/Ast_DCLK1/Ast_DCLK1.rds"

seurat_obj <- readRDS(in_path)
print(seurat_obj)
# metadata = seurat_obj@meta.data

# filter cells
filtered_cells <- rownames(seurat_obj@meta.data[!is.na(seurat_obj@meta.data$bmi_lv), ])
seurat_obj = subset(seurat_obj, cells = filtered_cells)
print(seurat_obj)

# normalize bmi (subtract mean, divide by SD)
bmi_lv = seurat_obj@meta.data$bmi_lv
mean = mean(bmi_lv)
sd=sd(bmi_lv)
print(paste0("Mean BMI: ", mean, ", SD BMI: ", sd))

bmi_normalized = (bmi_lv - mean) / sd
print(bmi_normalized)
seurat_obj@meta.data$bmi_normalized = bmi_normalized
head(seurat_obj@meta.data)

out_path = sub(".rds", "_bmi.rds", in_path)
saveRDS(seurat_obj, out_path)
print(paste0("Output saved to: ", out_path))
