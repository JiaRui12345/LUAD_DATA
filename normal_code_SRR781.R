rm(list = ls())
library(Seurat)##version 4.0.5
library(dplyr)
library(future)
library(future.apply)
library(DoubletFinder)
plan("multicore", workers = 18) ###compute cores
options(future.globals.maxSize = 50000 * 1024^2)

getwd()
#setwd("../../Data/10X/")
#list.files(path = "./",pattern = "^SRR")
####SRR781####
#SRR781<-Read10X(data.dir = "SRR81_outs/filtered_feature_bc_matrix/")
SRR781[1:5,1:5]

SRR781_object<- CreateSeuratObject(counts = SRR781, project = "normal", min.cells = 0, min.features = 200)
SRR781_object[["percent.mt"]] <- PercentageFeatureSet(SRR781_object, pattern = "^MT-")
hist(SRR781_object[["percent.mt"]]$percent.mt)
VlnPlot(SRR781_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(SRR781_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SRR781_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
SRR781_object
##select the cells with 300 genes at least and 4000 at most, the percent of mitochondrion genes is less than 10%
SRR781_val<- subset(SRR781_object, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)

SRR781_val@meta.data[1:5,]


#colname<-paste("SRR781_",colnames(SRR781_val),sep="")
#colname

#SRR781_val<-RenameCells(object = SRR781_val,colname)
#rm(list = c("BC2","bc2_object"))
SRR781_val@meta.data[1:5,]
SRR781_val$tissue_type<-"normal"


SRR781_val <- NormalizeData(SRR781_val, normalization.method = "LogNormalize", scale.factor = 10000)
SRR781_val <- FindVariableFeatures(SRR781_val, selection.method = "vst", nfeatures = 2000)

###scaling the data###
all.genes <- rownames(SRR781_val)
SRR781_val <- ScaleData(SRR781_val,vars.to.regress = c("percent.mt"))
SRR781_val
###perform linear dimensional reduction###
SRR781_val <- RunPCA(SRR781_val, features = VariableFeatures(object = SRR781_val))
print(SRR781_val[["pca"]], dims = 1:5, nfeatures = 5)
#dev.off()
VizDimLoadings(SRR781_val, dims = 1:2, reduction = "pca")


ElbowPlot(SRR781_val,ndims = 50)


####cluster the cells###
SRR781_val <- FindNeighbors(SRR781_val, dims = 1:30)
####Run non-linear dimensional reduction (UMAP/tSNE)###
SRR781_val <- RunUMAP(SRR781_val, dims = 1:30)
SRR781_val<-RunTSNE(SRR781_val,dims=1:30)
SRR781_val <- FindClusters(SRR781_val, resolution = 0.5)

DimPlot(SRR781_val, reduction = "umap",label = T)
#?DimPlot
DimPlot(SRR781_val, reduction = "tsne",label = T)




SRR781.markers <- FindAllMarkers(SRR781_val, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(SRR781.markers,"./markers.csv")

#marker <- FindMarkers(SRR781_val, ident.1 = "3",ident.2 = "7",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#肺癌注释基因
#epithelial cell,查一下这些上皮细胞所分布的群"CDH1","KRT18",
FeaturePlot(SRR780_val,features = c("EPCAM","KRT19","KRT7","KRT8"),reduction = "tsne",label = T)
#Fibroblasts
FeaturePlot(SRR780_val,features = c("DCN","THY1","COL1A1","COL1A2"),reduction = "tsne",label = T)
#Endothelial cells
FeaturePlot(SRR780_val,features = c("PECAM1","CLDN5","FLT1","RAMP2"),reduction = "tsne",label = T)
#T lymphocytes #"CD8A"
FeaturePlot(SRR780_val,features = c("CD3D","CD3E","CD3G","TRAC"),reduction = "tsne",label = T)
#NK cells
FeaturePlot(SRR780_val,features = c("NKG7","GNLY","NCAM1","KLRD1"),reduction = "tsne",label = T)
#B lymphocytes "IGHA2""IGHM"
FeaturePlot(SRR780_val,features = c("CD79A","IGHG3","MS4A1","CD22"),reduction = "tsne",label = T)
#Myeloid cells
FeaturePlot(SRR780_val,features = c("LYZ","MARCO","CD68","FCGR3A"),reduction = "tsne",label = T)
#MAST cells
FeaturePlot(SRR780_val,features = c("KIT","MS4A2","GATA2"),reduction = "tsne",label = T)


Idents(SRR781_val)<-SRR781_val$RNA_snn_res.0.5
SRR781_val<-RenameIdents(SRR781_val,"7"="Epithelial","16"="Epithelial","15"="Epithelial","13"="Epithelial",
                         "9"="Fibroblast","17"="Fibroblast",
                         "10"="MAST",
                         "18"="B","2"="T/NK","3"="T/NK","0"="T/NK",
                         "5"="Myeloid","8"="Myeloid","12"="Myeloid","14"="Myeloid","6"="Myeloid","4"="Myeloid","1"="Myeloid",
                         "19"="Endothelial","11"="Endothelial","20"="Endothelial","21"="Unknow","22"="Unknow")
DimPlot(SRR781_val,reduction = "tsne",label = F)
DimPlot(SRR781_val,reduction = "umap",label = F)
table(Idents(SRR781_val))###查看有多少个细胞

#SRR781_val$cell_cluster_origin<-Idents(SRR781_val)
SRR781_val@meta.data[1:5,]

#save.image(file="./SRR781_outs/SRR781_val.RData")
#load("./SRR781_outs/SRR781_val.RData")

SRR781_val$cluster<-Idents(SRR781_val)


#saveRDS(SRR781_val,file = "./SRR81_outs/SRR781_val.rds")
#SRR781_val<-readRDS("./SRR81_outs/SRR781_val.rds")

####Doublets identification#
library(DoubletFinder)
SRR781_val_db<-SRR781_val
sweep.res.list_SRR781_val <- paramSweep_v3(SRR781_val_db, PCs = 1:30)
sweep.stats_SRR781_val <- summarizeSweep(sweep.res.list_SRR781_val, GT = FALSE)
bcmvn_SRR781_val<-find.pK(sweep.stats_SRR781_val)
pK_value <- as.numeric(as.character(bcmvn_SRR781_val$pK[bcmvn_SRR781_val$BCmetric == max(bcmvn_SRR781_val$BCmetric)]))
annotations <- SRR781_val_db@meta.data$RNA_snn_res.0.5
homotypic.prop <- modelHomotypic(annotations)  ####get the estimated number
nExp_poi <- round(homotypic.prop*length(SRR781_val_db@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
SRR781_val_db <- doubletFinder_v3(SRR781_val_db, PCs = 1:30, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
SRR781_val_db <- doubletFinder(SRR781_val_db, PCs = 1:30, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
SRR781_val_db@meta.data[1:5,]

SRR781_val$doublets<-SRR781_val_db$DF.classifications_0.25_0.28_3758
SRR781_val@meta.data[1:5,]
setwd("F:/yufenxi/normal")
####saveRDS###
saveRDS(SRR781_val,file = "./SRR781_val.rds")
DimPlot(SRR781_val,reduction = "tsne",group.by = "doublets")
save.image(file="./normal.RData")


