library(nebula)
library(Seurat)
library(SeuratDisk)

indir <- "/home/anna_y/data/test/tiny_AD427_ADMR/"
filename <- "tiny_AD427_ADMR.rds"
argv <- commandArgs(T)
name <- "tiny_AD427_ADMR"
sample.col <- "Sample"
out_dir <- "/home/anna_y/data/test/deg/"

print("Loading data...")
print(paste0("input directory:", indir, filename))
print(paste0("output directory:", out_dir))

seurat_obj <- readRDS(paste0(indir, filename))
print(seurat_obj)
head(seurat_obj@meta.data)

deg.nebula <- function(Seurat_Obj, pathology, sample.col,
                       offset="total_counts",
                       ncore=12,
                       cpc=0, reml=1) {

  covariates = c("msex", "pmi", "LibType", "total_counts", "nFeature_RNA")
  covariates = covariates[covariates %in% colnames(Seurat_Obj@meta.data)]

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

  design <- model.matrix(as.formula(paste0("~", paste0(c(pathology, covariates), collapse="+"))),
                         data = Seurat_Obj@meta.data)

  neb <- tryCatch({
    nebula(seuratdata$count,
           seuratdata$id,
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

  neb$summary$FDR <- p.adjust(neb$summary[[paste0("p_", pathology)]], "fdr")
  neb$summary$log2FC <- neb$summary[[paste0("logFC_", pathology)]]

  neb_df <- neb$summary
  ovr <- neb$overdispersion
  colnames(ovr) <- paste0("overdispersion_", colnames(ovr))
  neb_df <- cbind(neb_df, ovr)
  neb_df$convergence <- neb$convergence
  neb_df$algorithm <- neb$algorithm
  rownames(neb_df) <- neb_df$gene

  neb_df <- neb_df[order(neb_df[[paste0("p_", pathology)]]),]
  write.table(neb_df, paste0(out_dir, name, ".", pathology, ".All.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  results <- neb_df[, c("gene", paste0("p_", pathology), paste0("logFC_", pathology), "FDR", "log2FC")]
  write.table(results, paste0(out_dir, name, ".", pathology, ".Clean.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

}

var <- "nCount_RNA"
deg.nebula(seurat_obj, var, sample.col)
print(paste("Running deg.nebula on", var))
