library(nebula)
library(Seurat)

argv <- commandArgs(T)
#name="Exc.DG.granule.cells"
name=argv[1]
pathology=argv[2]
indir=argv[3]
out_dir=argv[4]
sample.col="Sample"

print("Loading data")

Seurat_Obj<-readRDS(paste0(indir,"rna4.AD427_only.",name,".rds"))

#Seurat_Obj@meta.data[[pathology]] <- round(Seurat_Obj@meta.data[[pathology]], 1)
Seurat_Obj@meta.data[["pmi"]] <- round(Seurat_Obj@meta.data[["pmi"]], 0)

deg.nebula <- function(Seurat_Obj, pathology, sample.col,
                       offset="total_counts",
                       ncore=12) {

    covariates=c("msex","pmi","ADdiag3types","total_counts","nFeaturess_RNA")
    # covariates=c("msex","pmi","total_counts","nFeaturess_RNA")
    covariates = covariates[covariates %in%  colnames(Seurat_Obj@meta.data)]


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
    write.table(neb_df,paste0(out_dir,name,".",pathology,".All.tsv"),
                sep = "\t", row.names = FALSE, quote = F)
    results=neb_df[,c("gene",paste0("p_", pathology),paste0("logFC_", pathology),"FDR","log2FC")]
    write.table(results,paste0(out_dir,name,".",pathology,".Clean.tsv"),
                sep = "\t", row.names = FALSE, quote = F)
}


deg.nebula(Seurat_Obj,pathology,sample.col)
