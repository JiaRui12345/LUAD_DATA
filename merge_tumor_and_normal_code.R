####全体样本的整合分析####
####merge the scRNA-data####
setwd("F:/yufenxi/merge_tumor_and_normal")
rm(list=ls())
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(dplyr)
library(msigdbr)
library(clusterProfiler)
#BiocManager::install("clusterProfiler")
#BiocManager::install("msigdbr")
#设置电脑线程和内存
plan("multicore", workers = 18) ###set the compute core
options(future.globals.maxSize = 50000 * 1024^2)
getwd()
#读入数据 SRR781为正常，SRR780为肿瘤
#setwd("~/Documents_PC/scRNA-seq/Data/10X")
#SRR780<-readRDS("./SRR80_outs/SRR780_val.rds")
#SRR781<-readRDS("./SRR81_outs/SRR781_val.rds")
#SRR782<-readRDS("./SRR82_outs/SRR782_val.rds")
#SRR783<-readRDS("./SRR83_outs/SRR783_val.rds")
#查看一下每一个样本有多少个细胞
SRR780@meta.data[1:5,]
SRR781@meta.data[1:5,]
#SRR782@meta.data[1:5,]
#SRR783@meta.data[1:5,]
###########
#做成一个大的seurat文件
sce.mergeTEN<- merge(SRR781,y=c(SRR780),project = "scTEN")
#查看一下，一定要有4个样本的normal这样的分组信息
sce.mergeTEN@meta.data[1:5,]
rm(SRR780,SRR781)
gc()
saveRDS(sce.mergeTEN,file="./sce.mergeTEN_2.rds")
#sce.mergeTEN<-readRDS(file="./sce.mergeTEN.rds")

####integration with harmony####
library(devtools)
#install_github("immunogenomics/harmony")
library(harmony)
gc()
sce.mergeTEN@meta.data[1:5,]
#细胞不同增殖状态可能会影响分组，为了避免这种影响而进行这一步分析
#s期的基因
s.genes <- cc.genes$s.genes
#g2m期的基因
g2m.genes <- cc.genes$g2m.genes
#推算这个细胞所处的细胞分期
sce.mergeTEN <- CellCycleScoring(sce.mergeTEN, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
sce.mergeTEN@meta.data[1:5,]
#进行标准化
sce.mergeTEN<-NormalizeData(sce.mergeTEN,verbose = T) 
#选择2000个top基因进行后续分析
sce.mergeTEN<-FindVariableFeatures(sce.mergeTEN,selection.method = "vst", nfeatures = 2000)
#去除线粒体,S期,G2M期的影响
sce.mergeTEN<-ScaleData(sce.mergeTEN,vars.to.regress = c("percent.mt","S.Score","G2M.Score"),verbose = FALSE)
#pca降维
sce.mergeTEN<-RunPCA(sce.mergeTEN,verbose = T,npcs = 50)
ElbowPlot(sce.mergeTEN,ndims = 50)#有批次效应时，批次效应选的尽量大
#查看有无批次效应
p1 <- DimPlot(object = sce.mergeTEN, reduction = "pca", pt.size = .1, group.by = "orig.ident",raster=FALSE)
p2 <- VlnPlot(object = sce.mergeTEN, features = "PC_1", group.by = "orig.ident", pt.size = .1,raster=FALSE)
CombinePlots(plots=list(p1,p2))
#去除因为来自不同个体的批次效应
sce.mergeTEN<-RunHarmony(sce.mergeTEN,group.by.vars = c("orig.ident"), plot_convergence = TRUE)
harmony_embeddings <- Embeddings(sce.mergeTEN, 'harmony')
dim(harmony_embeddings)
#再次查看处理后的批次效应情况
p3 <- DimPlot(object = sce.mergeTEN, reduction = "harmony", pt.size = .1, group.by = "orig.ident",raster=FALSE)
p4 <- VlnPlot(object = sce.mergeTEN, features = "harmony_1", group.by = "orig.ident", pt.size = .1,raster=FALSE)
CombinePlots(plots=list(p3,p4))
p1+p3
#跑ump和tsne
sce.mergeTEN <- sce.mergeTEN %>% 
  RunUMAP(reduction = "harmony", dims = 1:50) %>% 
  RunTSNE(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50)
#进行分群
sce.mergeTEN<-FindClusters(sce.mergeTEN,resolution = 0.8)#大样本想要分得更细，resolution可以调的大一点
table(Idents(sce.mergeTEN))
#Idents(sce_inter)<-sce_inter$seurat_clusters
sce.mergeTEN@meta.data[1:5,]
#meta.data <- sce.mergeTEN@meta.data

#画分群的
p1<-DimPlot(sce.mergeTEN,reduction = "tsne",label = T,raster=FALSE)
p1
DimPlot(sce.mergeTEN,reduction = "umap",label = T,raster=FALSE)
#DimPlot(sce.mergeTEN,reduction = "tsne",label = T,split.by = "tissue_type",raster=FALSE)
#DimPlot(sce.mergeTEN,reduction = "umap",label = T,split.by = "tissue_type",raster=FALSE)
table(sce.mergeTEN@meta.data$old.ident)
#查看前期分群在本次分群中的效果
sce.mergeTEN1@meta.data[1:5,]
p2<-DimPlot(sce.mergeTEN,reduction = "tsne",group.by = "old.ident",label = F,raster=FALSE)
p2
DimPlot(sce.mergeTEN,reduction = "tsne",group.by = "old.ident",split.by = "tissue_type",label = F,raster=FALSE)
DimPlot(sce.mergeTEN,reduction = "tsne",group.by = "old.ident",split.by = "doublets",label = F,raster=FALSE)

#DimPlot(sce.mergeTEN,group.by = "cell_cluster_origin",reduction = "tsne",label = T)
table(Idents(sce.mergeTEN),sce.mergeTEN$old.ident)

DimPlot(sce.mergeTEN,reduction = "tsne",label = T)

####通过经典marker进行标注####
#肺癌注释基因
#epithelial cell,查一下这些上皮细胞所分布的群"CDH1","KRT7",
FeaturePlot(sce.mergeTEN,features = c("EPCAM","CDH1","KRT7","KRT8"),reduction = "tsne",label = T,raster=FALSE)
#Fibroblasts
FeaturePlot(sce.mergeTEN,features = c("DCN","THY1","COL1A1","COL1A2"),reduction = "tsne",label = T,raster=FALSE)
#Endothelial cells
FeaturePlot(sce.mergeTEN,features = c("PECAM1","CLDN5","FLT1","RAMP2"),reduction = "tsne",label = T,raster=FALSE)
#T lymphocytes"TRAC","CD8A",
FeaturePlot(sce.mergeTEN,features = c("CD3D","CD3E","CD3G","TRAC"),reduction = "tsne",label = T,raster=FALSE)
#NK cells
FeaturePlot(sce.mergeTEN,features = c("NKG7","GNLY","NCAM1","KLRD1"),reduction = "tsne",label = T,raster=FALSE)
#B lymphocytes"MS4A1","CD22""IGHG3","IGHA2"
FeaturePlot(sce.mergeTEN,features = c("CD79A","IGHM","MS4A1","CD22"),reduction = "tsne",label = T,raster=FALSE)
#Myeloid cells
FeaturePlot(sce.mergeTEN,features = c("LYZ","MARCO","CD68","FCGR3A"),reduction = "tsne",label = T,raster=FALSE)
#MAST cells
FeaturePlot(sce.mergeTEN,features = c("KIT","MS4A2","GATA2"),reduction = "tsne",label = T,raster=FALSE)

#进行findallmarkers分析
#sce.mergeTEN.markers <- FindAllMarkers(sce.mergeTEN, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#write.csv(sce.mergeTEN.markers,"../10X/sce.mergeTEN.markers.csv")

#marker <- FindMarkers(sce.mergeTEN, ident.1 = "21",ident.2 = "8",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#c<-table(Idents(sce.mergeTEN),sce.mergeTEN$cell_predict)
#c
#write.csv(c,file="../10X/cc.csv")
#p1+p2
#save.image(file="./all.RData")
#进行亚群标注
sce.mergeTEN<-RenameIdents(sce.mergeTEN,"5"="Epithelial","16"="Epithelial","19"="Epithelial","21"="Epithelial","23"="Epithelial",
                           "24"="Epithelial","25"="Epithelial","26"="Epithelial","27"="Epithelial","29"="Epithelial","33"="Epithelial",
                           "0"="T/NK","1"="T/NK","2"="T/NK","9"="T/NK","18"="T/NK","22"="T/NK",
                           "4"="B","15"="B","3"="Myeloid","6"="Myeloid","7"="Myeloid","8"="Myeloid","11"="Myeloid",
                           "13"="Myeloid","17"="Myeloid","20"="Myeloid","12"="MAST",
                           "10"="Fibroblast","30"="Fibroblast",
                           "14"="Endothelial","32"="Endothelial",
                           "31"="Unknown","29"="Unknown","28"="Unknown")
#画图
DimPlot(sce.mergeTEN,reduction = "tsne",label = F,raster=FALSE)
DimPlot(sce.mergeTEN,reduction = "umap",label = F,raster=FALSE)
#table(Idents(SC294_val))
DimPlot(sce.mergeTEN,reduction = "tsne",label = F,group.by = "doublets",raster=FALSE)
DimPlot(sce.mergeTEN,reduction = "tsne",label = F,split.by = "orig.ident",raster=FALSE)

sce.mergeTEN$Seurat_harmony<-Idents(sce.mergeTEN)
sce.mergeTEN@meta.data[1:5,]
#安装上述刚标记的cluster进行findallmarkers分析
#markers<-FindAllMarkers(sce.mergeTEN,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
#setwd("G:/yufenxi/all")
#write.csv(markers,file="./markers.csv")#这个表可以放在附件里
save.image(file="./merge.RData")


#####对图进行美化####
#BiocManager::install("ggsci")
library(ggsci)
##define the color
??ggsci
#寻找一个合适的配色模板
#pal_npg()(38)
cors<-pal_igv()(38)
cors<- pal_aaas()(38)
#cors <- c(cors,"grey")
features = c("EPCAM","CDH1","KRT8","KRT7",
             "DCN","COL1A2","COL1A1","THY1",
             "CLDN5","RAMP2","PECAM1","FLT1",
             "CD3D","CD3E","CD3G","TRAC","NKG7","GNLY",
             "CD79A","IGHG3","MS4A1",
             "CD68","LYZ","FCGR3A","MARCO",
             "MS4A2","GATA2","KIT")

length(features)
VlnPlot(sce.mergeTEN,features = features[1:8],cols = cors,pt.size = 0,raster=FALSE,ncol = 4)
VlnPlot(sce.mergeTEN,features = features[9:16],cols = cors,pt.size = 0,raster=FALSE,ncol = 4)
VlnPlot(sce.mergeTEN,features = features[17:24],cols = cors,pt.size = 0,raster=FALSE,ncol = 4)
VlnPlot(sce.mergeTEN,features = features[25:28],cols = cors,pt.size = 0,raster=FALSE,ncol = 4)

VlnPlot(sce.mergeTEN,features = c("CD14","VEGFA"),cols = cors,pt.size = 0)

DotPlot(sce.mergeTEN,features = features,cols = c("grey","darkblue"))

#####画一个marker基因平均表达的heatmap图#####
#cors
####matrix
genemeanMatrix<-AverageExpression(sce.mergeTEN)
genemeanMatrix<-genemeanMatrix$RNA
#感兴趣的marker基因
features_to_anno = c("EPCAM","CDH1","KRT8","KRT7",
                     "DCN","COL1A2","COL1A1","THY1",
                     "CLDN5","RAMP2","PECAM1","FLT1",
                     "CD3D","CD3E","CD3G","TRAC","NKG7","GNLY",
                     "CD79A","IGHG3","MS4A1",
                     "CD68","LYZ","FCGR3A","MARCO",
                     "MS4A2","GATA2","KIT")
length(features_to_anno)
geneMatrix<-genemeanMatrix[features_to_anno,]
#rownames
cc<-rownames(geneMatrix)
cc
#cc[which(!(cc%in%features_to_anno))]<-""
#cc
table(Idents(sce.mergeTEN))
levels(Idents(sce.mergeTEN))##知道细胞亚群包含哪些
#information for columns（根据之前的标记）
annotation_col = data.frame(
  CellType = c("Epithelial", "T/NK","B","Myeloid","MAST","Fibroblast","Endothelial", "Unknown")
)
rownames(annotation_col) = colnames(geneMatrix)
#define the seq
annotation_col$CellType<-factor(annotation_col$CellType,levels= c("Epithelial", "T/NK","B","Myeloid","MAST","Fibroblast","Endothelial", "Unknown"))
levels(Idents(sce.mergeTEN))
cors
paste(c("Epithelial", "T/NK","B","Myeloid","MAST","Fibroblast","Endothelial","Unknown"),cors,sep='"="')
#color to use in heatmap
ann_colors = list(
  CellType = c("Epithelial"="#5050FFFF", "T/NK"="#CE3D32FF",        "B"="#749B58FF",    
               "Myeloid"="#466983FF",     "MAST"="#BA6338FF",        "Fibroblast"="#5DB1DDFF",  "Endothelial"="#802268FF",
               "Unknown"="#6BD76BFF")
  #GeneClass = c("CD4+" = "#7570B3", "CD8+" = "#E7298A","NK" = "#66A61E")
)
?pheatmap
library(RColorBrewer)
#dev.off()
#dev.new()
ComplexHeatmap::pheatmap(as.matrix(geneMatrix), annotation_col = annotation_col,cluster_rows = F,labels_row  = cc,
                         annotation_colors = ann_colors,cluster_cols = F,scale = "row",color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(50))
#gaps_row = c(8, 14) )

