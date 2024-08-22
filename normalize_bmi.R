library(Seurat)

filter_cells = function(seurat_obj, var) {
    # return a subset of cells that contain the given variable
    filtered_cells <- rownames(seurat_obj@meta.data[!is.na(seurat_obj@meta.data[[var]]), ])
    seurat_obj = subset(seurat_obj, cells = filtered_cells)
    return (seurat_obj)
}

# Function to standardize BMI values in the Seurat object
standardize_bmi = function(seurat_obj, save=NULL) {
    # Extract the BMI values
    bmi_lv = seurat_obj@meta.data$bmi_lv

    # Check for non-missing BMI values
    non_missing_bmi = !is.na(bmi_lv)

    # Calculate mean and SD only for non-missing BMI values
    mean_bmi = mean(bmi_lv[non_missing_bmi])
    sd_bmi = sd(bmi_lv[non_missing_bmi])
    print(paste0("Mean BMI: ", mean_bmi, ", SD BMI: ", sd_bmi))

    # Normalize BMI (subtract mean, divide by SD) only for non-missing values
    bmi_normalized = bmi_lv
    bmi_normalized[non_missing_bmi] = (bmi_lv[non_missing_bmi] - mean_bmi) / sd_bmi

    # Store the normalized BMI in the metadata
    seurat_obj@meta.data$bmi_normalized = bmi_normalized
    print(head(seurat_obj@meta.data))

    # Save the modified Seurat object if a path is provided
    if (!is.null(save)) {
        saveRDS(seurat_obj, save)
        print(paste0("Output saved to: ", save))
    }

    return(seurat_obj)
}

if (sys.nframe() == 0) {

    in_path = commandArgs(T)[1] # path to .rds file
    seurat_obj <- readRDS(in_path)
    seurat_obj <- standardize_bmi(seurat_obj, save=in_path) # overwrite the input .rds file

}
