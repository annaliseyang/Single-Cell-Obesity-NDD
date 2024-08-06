library(nebula)
library(Seurat)
library(SeuratDisk)

# indir="/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_ADMR_multiome_snRNA_snATAC/01_multiome_snRNA/03_DEG_Nebula/ClassBrainRegion/"
indir = "/home/anna_y/data/write/"
outdir = "/home/anna_y/data/write/"
filename = "AD427_ADMR_meta_Jul22_2024.h5seurat"
Convert(indir, dest = "h5seurat", overwrite = TRUE)

argv <- commandArgs(T)
#name="Exc.DG.granule.cells"
name=argv[1]
sample.col="Sample"
out_dir="/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_ADMR_multiome_snRNA_snATAC/01_multiome_snRNA/03_DEG_Nebula/ClassBrainRegion_gpath/"

print("Loading data")

# Seurat_Obj<-readRDS(paste0(indir,"rna9.ADMR.All.",name,".rds"))
# Seurat_Obj <- readRDS(paste0(indir, filename))
Seurat_Obj <- LoadH5Seurat(paste0(indir, filename))

deg.nebula <- function(Seurat_Obj, pathology, sample.col, 
                       offset="total_counts", 
                       ncore=12,
                       cpc=0, reml=1) {
  
  covariates=c("msex","pmi","LibType","total_counts","nFeature_RNA")
  covariates = covariates[covariates %in%  colnames(Seurat_Obj@meta.data)]
  
  #pathology="Epigenetic_Information"
  
  seuratdata <- scToNeb(obj = Seurat_Obj, 
                        assay = "RNA", 
                        id = sample.col, 
                        pred = c(covariates,pathology), 
                        offset="total_counts")
  
  design = model.matrix(as.formula(paste0("~", paste0(c(pathology, covariates), collapse="+"))), 
                        data=Seurat_Obj@meta.data)
  
  neb = nebula(seuratdata$count,
               seuratdata$id,
               pred=design,
               model="NBGMM", 
               ncore=ncore,
               offset=seuratdata$offset)
  
  ### Now to parse output
  neb$summary$FDR = p.adjust(neb$summary[[paste0("p_", pathology)]], "fdr")
  neb$summary$log2FC = neb$summary[[paste0("logFC_", pathology)]]
  
  #diff=neb$summary[neb$summary$p_Epigenetic_Information<0.01,]
  
  neb_df = neb$summary
  ovr = neb$overdispersion
  colnames(ovr) = paste0("overdispersion_", colnames(ovr))
  #print(str(ovr))
  neb_df = cbind(neb_df, ovr)
  neb_df$convergence = neb$convergence
  neb_df$algorithm = neb$algorithm
  rownames(neb_df) = neb_df$gene
  
  neb_df<-neb_df[order(neb_df[[paste0("p_", pathology)]]),]
  write.table(neb_df,paste0(out_dir,name,".",pathology,".All.tsv"), sep = "\t", row.names = F,quote = F)
  results=neb_df[,c("gene",paste0("p_", pathology),paste0("logFC_", pathology),"FDR","log2FC")]
  write.table(results,paste0(out_dir,name,".",pathology,".Clean.tsv"), sep = "\t", row.names = F,quote = F)
  
}


deg.nebula(Seurat_Obj,"gpath",sample.col)