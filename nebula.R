library(nebula)
library(Seurat)
library(SeuratDisk)

var <- "bmi_lv"
# var <- "AD_states"
indir = commandArgs(T)[1] # eg. /home/anna_y/data/write/Class/Ast/

filename <- list.files(indir, pattern=".rds")[1]

while (is.na(filename) || length(filename) == 0) {
  print("No .rds file found in the input directory. Checking again in 5 minutes.")
  Sys.sleep(300) # wait 5 minutes before checking again
  filename <- list.files(indir, pattern=".rds")[1]
}

# Proceed with the rest of the code after an .rds file is found
print(paste("Found .rds file:", filename))

name <- sub(".rds", "", filename)
sample.col <- "Sample"
# out_dir <- "~/results/deg/"
out_dir <- sub("data/write", paste0("results/", "deg_", var), indir)

# create the out_dir if it doesn't exist
dir.create(out_dir, recursive = TRUE)

print("Loading data...")
print(paste0("input directory: ", file.path(indir, filename)))
print(paste0("output directory: ", out_dir))

seurat_obj <- readRDS(file.path(indir, filename))
print(seurat_obj)
# head(seurat_obj@meta.data)
print("Metadata dimensions:")
print(dim(seurat_obj@meta.data))

deg.nebula <- function(Seurat_Obj, pathology, sample.col,
                       offset = "total_counts",
                       ncore = 12,
                       cpc = 0, reml = 1) {

  covariates = c("msex", "pmi", "total_counts", "nFeature_RNA", "age_death", "gpath")
  covariates = covariates[covariates %in% colnames(Seurat_Obj@meta.data)]
  # Seurat_Obj@meta.data$pathology <- as.factor(Seurat_Obj@meta.data$pathology) # convert to factor
  head(Seurat_Obj@meta.data$pathology)

  print("Covariates:")
  print(covariates)

  # filtered_meta_data <- Seurat_Obj@meta.data[!is.na(Seurat_Obj@meta.data[[pathology]]), ]
  # Seurat_Obj <- subset(Seurat_Obj, cells = rownames(filtered_meta_data))
  metadata <- Seurat_Obj@meta.data

  print("Filtering cells...")
  filtered_cells <- rownames(metadata[!is.na(metadata$bmi_lv), ])
  print(paste0("Number of cells after filtering by ", pathology, ": ", length(filtered_cells)))
  Seurat_Obj <- subset(Seurat_Obj, cells = filtered_cells)
  print(Seurat_Obj)
  metadata <- metadata[rownames(Seurat_Obj), ]

  # Convert Seurat object to Nebula data format
  seuratdata <- tryCatch({
    scToNeb(obj = Seurat_Obj,
            assay = "RNA",
            id = sample.col,
            pred = c(covariates, pathology),
            offset = "total_counts")
  }, error = function(e) {
    cat("Error in scToNeb: ", e$message, "\n")
    NULL
  })

  if (is.null(seuratdata)) {
    stop("Failed to convert Seurat object to Nebula data. Exiting...")
  }

  # print("Grouping data...")
  # # Group cells by sample ID
  # data_g <- group_cell(count = seuratdata$count, id = seuratdata$id, pred = seuratdata$pred)
  # print("Data grouped by subject id!")
  data_g <- seuratdata
  print(dim(data_g$count))

  # Convert predictor matrix to dataframe
  df <- as.data.frame(data_g$pred)
  # df <- as.data.frame(metadata)
  colnames(df) <- c(covariates, pathology)

  print("Creating formula...")
  formula <- as.formula(paste0("~", paste0(c(pathology, covariates), collapse = "+")))
  print(formula)
  design <- model.matrix(formula, data = df)

  print("Number of rows in design matrix:")
  print(nrow(design))
  print("Number of columns in grouped count matrix:")
  print(ncol(data_g$count))

  # Find common cells between design and count matrices
  common_cells <- intersect(rownames(df), colnames(data_g$count))
  print("Number of common cells:")
  print(length(common_cells))

  if (length(common_cells) == 0) {
    stop("No common cells found between design matrix and count matrix.")
  }

  # Subset both design matrix and count matrix to common cells
  design <- design[common_cells, , drop = FALSE]
  data_g$count <- data_g$count[, common_cells, drop = FALSE]

  print("Number of rows in design matrix after subsetting:")
  print(nrow(design))
  print("Number of columns in count matrix after subsetting:")
  print(ncol(data_g$count))

  if (nrow(design) == 0 || ncol(data_g$count) == 0) {
    stop("Error: Design matrix or count matrix is empty after filtering. Please check your data.")
  }


  # Run nebula
  neb <- tryCatch({
    nebula(data_g$count,
           data_g$id,
           pred = design,
           model = "NBGMM",
           ncore = ncore,
           offset = seuratdata$offset)
  }, error = function(e) {
    cat("Error in nebula: ", e$message, "\n")
    NULL
  })

  if (is.null(neb)) {
    stop("Failed to run nebula. Exiting...")
  }

  # Process and save results
  neb$summary$FDR <- p.adjust(neb$summary[[paste0("p_", pathology)]], "fdr")
  neb$summary$log2FC <- neb$summary[[paste0("logFC_", pathology)]]

  neb_df <- neb$summary
  ovr <- neb$overdispersion
  colnames(ovr) <- paste0("overdispersion_", colnames(ovr))
  neb_df <- cbind(ovr, neb_df)
  neb_df$convergence <- neb$convergence
  neb_df$algorithm <- neb$algorithm
  rownames(neb_df) <- neb_df$gene

  neb_df <- neb_df[order(neb_df[[paste0("p_", pathology)]]), ]

  out_path_all <- file.path(out_dir, paste0(name, ".", pathology, ".All.tsv"))
  write.table(neb_df, out_path_all, sep = "\t", row.names = FALSE, quote = FALSE)
  print(paste("All results saved in:", out_path_all))

  results <- neb_df[, c("gene", paste0("p_", pathology), paste0("logFC_", pathology), "FDR", "log2FC")]
  out_path_clean <- file.path(out_dir, paste0(name, ".", pathology, ".Clean.tsv"))
  write.table(results, out_path_clean, sep = "\t", row.names = FALSE, quote = FALSE)
  print(paste("Clean results saved in:", out_path_clean))
}


print(paste("Running deg.nebula on", var))
deg.nebula(seurat_obj, var, sample.col)
