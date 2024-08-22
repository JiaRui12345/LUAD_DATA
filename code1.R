library(Seurat)	
library(dyno)
library(monocle3)	
library(ggplot2)	
library(dplyr)	
library(CellChat)	
#devtools::install_github("sqjin/CellChat")
#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'S4Vectors', 'SingleCellExperiment','SummarizedExperiment', 'batchelor', 'Matrix.utils'))
#install.packages("devtools")
#library(speedglm)
#setwd("C:/Program Files/R/R-4.2.3/library/speedglm_0.3-4")
#devtools::load_all("C:/Program Files/R/R-4.2.3/library/speedglm")
#devtools::install_version('speedglm', '0.3-4', repos = 'https://packagemanager.rstudio.com/cran/2023-03-31')
#devtools::install_github('cole-trapnell-lab/leidenbase',force = TRUE)
#devtools::install_github('cole-trapnell-lab/monocle3')
#devtools::install_github("dynverse/dyno")
seurat.obj <- readRDS("sce.mergeTEN.rds")	
seurat.obj <- sce.mergeTEN
#subcluster	
fibroblast <- subset(seurat.obj,Seurat_harmony=="Fibroblast")	
fibroblast <- FindVariableFeatures(fibroblast, selection.method = "vst", nfeatures = 2000)	
fibroblast <- RunPCA(fibroblast, features = VariableFeatures(object = fibroblast))	
fibroblast <- FindNeighbors(fibroblast, dims = 1:10)	
fibroblast <- FindClusters(fibroblast, resolution = 0.5)	
fibroblast <- RunUMAP(fibroblast, dims = 1:10)	
DimPlot(fibroblast, reduction = "umap")	
DimPlot(fibroblast, reduction = "umap",split.by = "tissue_type",label = F,raster=FALSE)	

#0  "MFAP5","PCOLCE2","PI16"
FeaturePlot(fibroblast,features = c("MFAP5","PCOLCE2","PI16"),reduction = "umap",label = T)
#1,6  "BMP5","LTBP2","TCF21"
FeaturePlot(fibroblast,features = c("BMP5","LBH","NPNT","TCF21"),reduction = "umap",label = T)
#2,5  "MYH11","PCSK1N","WIF1" COLLAGEN
FeaturePlot(fibroblast,features = c("COL11A1","COL10A1","COL1A2","FN1"),reduction = "umap",label = T)
#2
#FeaturePlot(fibroblast,features = c("FN1","COL1A2","TPM1","SPARC"),reduction = "umap",label = T)
#3  "IGHG4","IGLC3","IGHG1","IGHG3" ICAF
FeaturePlot(fibroblast,features = c("IGHG4","IGHG1","IGHG3","IGLC3"),reduction = "umap",label = F)
#4,8  "MYH11","ACTA2","TINAGL1","ACTG2","DES" MYCAF
FeaturePlot(fibroblast,features = c("TINAGL1","MYH11","ACTG2","TAGLN"),reduction = "umap",label = T)
#5   "POSTN","COL11A1","MMP11"
FeaturePlot(fibroblast,features = c("POSTN","MMP11","CTHRC1"),reduction = "umap",label = F)
#6  "FGFR4","FIGF","TCF21" 
#FeaturePlot(fibroblast,features = c("FGFR4","FIGF","RGCC","LIMCH1","TCF21"),reduction = "umap",label = T)
#7  "MYC","ADAMTS1","IL6" INCAT
FeaturePlot(fibroblast,features = c("MYC","IL6","CXCL2","NAMPT"),reduction = "umap",label = T)
#8  "MYH11","PCSK1N","WIF1"
#FeaturePlot(fibroblast,features = c("WIF1","FAM150A","PCSK1N","SCX","MYH11"),reduction = "umap",label = T)
#2,5  "MYH11","PCSK1N","WIF1"
#FeaturePlot(fibroblast,features = c("COL11A1","COL8A1","COL10A1","CTHRC1"),reduction = "umap",label = T)
#markers	
n_clust <- names(table(fibroblast@active.ident)[table(fibroblast@active.ident)!=0])	
marker.list <- list()	
for(i in n_clust){	
  ident1 <- i	
  markers <- FindMarkers(fibroblast,	
                         ident.1 = ident1, 	
                         only.pos = FALSE, 	
                         min.pct = 0.25, 	
                         logfc.threshold = 0.25)	
  marker.list[[i]] <- rownames(markers[1:20,])  	
}	

new_cluster <- c("COL14A1+","COL13A1+","COL1A2+",	
                 "Immune_related CAF","Myofibroblast","MyCAF",	
                 "COL13A1+","Inflammation_related","Myofibroblast")

names(new_cluster) <- levels(fibroblast)	
fibroblast <- RenameIdents(fibroblast,new_cluster)	
DimPlot(fibroblast, 	
        repel = TRUE,
        reduction = "umap", 
        label=TRUE, 
        pt.size = .1) + NoLegend()


saveRDS(fibroblast, file="fb_seurat.rds")