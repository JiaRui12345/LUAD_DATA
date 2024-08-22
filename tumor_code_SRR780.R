####首先认识数据结构####
###10XGenomics 的filtered_feature_bc_matrix
##https://www.ncbi.nlm.n[ih.gov/geo/query/acc.cgi?acc=GSE134520
###BD RSEC_MolsPerCell.csv
###SMART-seq RPKM.txt
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
install.packages("DoubletFinder")
rm(list = ls())

#install.packages("Seurat")
####首先我们先处理10X Genomics的数据
library(Seurat)##version 4.0.5
library(dplyr)
library(future)
library(future.apply)
library(DoubletFinder)
plan("multicore", workers = 18) ###compute cores，设置跑的线程
options(future.globals.maxSize = 50000 * 1024^2)
#getwd()
setwd("F:/yufenxi/tumor/")
#setwd("./10X/")

#setwd("~/Documents_PC/scRNA-seq/Data/10X")
#list.files(path = "./",pattern = "^SRR")

####数据导入和准备####
#SRR780<-Read10X(data.dir = "SRR80_outs/filtered_feature_bc_matrix/")
SRR780[1:5,1:5]

####数据处理和质量控制####
#将matrix整理成large seurat格式
SRR780_object<- CreateSeuratObject(counts = SRR780, project = "tumor", min.cells = 0, min.features = 200)
#计算并查看线粒体基因的占比
SRR780_object[["percent.mt"]] <- PercentageFeatureSet(SRR780_object, pattern = "^MT-")
hist(SRR780_object[["percent.mt"]]$percent.mt)
#nFeature_RNA表达数，nCount_RNA细胞数，percent.mt线粒体比例三者间的关系
VlnPlot(SRR780_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SRR780_object
plot1 <- FeatureScatter(SRR780_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SRR780_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
SRR780_object
##select the cells with 300 genes at least and 4000 at most, the percent of mitochondrion genes is less than 10%
###根据图1和图2 去除表达小于200或大于7500或>20%线粒体基因
SRR780_val<- subset(SRR780_object, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
#数据查看
SRR780_val@meta.data[1:5,]
dim(SRR780_val@meta.data)
colnames(SRR780_val)
rownames(SRR780_val@meta.data)
#若样本本身没有加入分组信息，对样本名字进行补全
#colname<-paste("SRR780_",colnames(SRR780_val),sep="")
#colname
#SRR780_val<-RenameCells(object = SRR780_val,colname)
#rm(list = c("BC2","bc2_object"))
SRR780_val@meta.data[1:5,]
#在meta_data里增加一列标记所有样本所属类型
SRR780_val$tissue_type<-"tumor"
SRR780_val@meta.data[1:5,]
#对样本进行标准化
SRR780_val <- NormalizeData(SRR780_val, normalization.method = "LogNormalize", scale.factor = 10000)
#使用样本中vst计算中表达差异最大的2000基因进行分群
SRR780_val <- FindVariableFeatures(SRR780_val, selection.method = "vst", nfeatures = 2000)
#查看前2000个基因是哪些和all.gene
SRR780_val@assays$RNA@var.features
length(SRR780_val@assays$RNA@var.features)

all.genes <- rownames(SRR780_val)
all.genes[1000:10010]
#SRR780_val <- ScaleData(SRR780_val, features = all.genes)
#去除线粒体基因的占比并进行归一化
SRR780_val <- ScaleData(SRR780_val,vars.to.regress = c("percent.mt"))
SRR780_val

####进行pca分析####
#对高突变的2000个基因进行pca降维
SRR780_val <- RunPCA(SRR780_val, features = VariableFeatures(object = SRR780_val))
print(SRR780_val[["pca"]], dims = 1:5, nfeatures = 5)
#dev.off()
#主成分中高变化基因（可以较好区分该成分发基因）从大到小的列表图
VizDimLoadings(SRR780_val, dims = 1:5, reduction = "pca")
#选择哪个主成分（看看选择哪个维度主成分时误差下降基于平稳）进行细胞分群
ElbowPlot(SRR780_val,ndims = 50)


####进行细胞聚类####
#根据上一步中的合适维度对细胞进行分群，一下三行代码中dims都是根据上一步进行调整，一般在20-30都差不多
SRR780_val <- FindNeighbors(SRR780_val, dims = 1:30)
####Run non-linear dimensional reduction (UMAP/tSNE)###
#umap降维
SRR780_val <- RunUMAP(SRR780_val, dims = 1:30)
#tsne降维
SRR780_val<-RunTSNE(SRR780_val,dims=1:30)
#进行分群，一般resolution为0.5-1，值越高亚群越多越细但也越有可能出现过拟合问题
SRR780_val <- FindClusters(SRR780_val, resolution = 0.5)
#画umap图，不要标记数字的话label取F
DimPlot(SRR780_val, reduction = "umap",label = T)
#画tsne图
DimPlot(SRR780_val, reduction = "tsne",label = T)
#查看各个分群细胞数量
table(Idents(SRR780_val))
#getwd()
#对所有亚群进行比较，查找各群中高表达，比例至少是25%，跟其他群比较logFC大于0.25（阈值越大越少）的基因
SRR780.markers <- FindAllMarkers(SRR780_val, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(SRR780.markers,"./markers.csv")
#也可以单个群之间进行比较，下方为比较3比7群获得表达高的基因
marker <- FindMarkers(SRR780_val, ident.1 = "3",ident.2 = "7",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


####进行人工注释的方法####
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

#根据上述分型进行标注相应亚型的类型
SRR780_val<-RenameIdents(SRR780_val,"21"="Epithelial","18"="Epithelial","26"="Epithelial","19"="Epithelial",
                         "24"="Epithelial","12"="Epithelial","11"="Epithelial","15"="Epithelial","4"="Epithelial",
                         "25"="Epithelial","8"="Fibroblast","16"="Endothelial","7"="MAST","10"="B","3"="B",
                         "20"="T/NK","0"="T/NK","1"="T/NK","6"="T/NK","9"="T/NK","22"="T/NK","17"="T/NK",
                         "2"="Myeloid","13"="Myeloid","14"="Myeloid","23"="Myeloid","5"="Myeloid","27"="Unknown","28"="Unknown")

DimPlot(SRR780_val,reduction = "tsne",label = F)
DimPlot(SRR780_val,reduction = "umap",label = F)
table(Idents(SRR780_val))###查看有多少个细胞
#saveRDS(SRR780_val,file = "./SRR780_val.rds")
#正常查看
SRR780_val@meta.data[1:5,]
table(Idents(SRR780_val))
#Idents(SRR780_val)<-SRR780_val$RNA_snn_res.0.5
#table(Idents(SRR780_val))
#对增加一列，把建立好的单细胞分群名称信息给到cluster列
#SRR780_val$cell_cluster_origin<-Idents(SRR780_val)
SRR780_val$cluster<-Idents(SRR780_val)

SRR780_val@meta.data[1:20,]
#saveRDS(SRR780_object,file = "./SRR780_obj.rds")
#save.image(file="./SRR780_dou.RData")
#load("./SRR780_outs/SRR780_val.RData")

table(Idents(SRR780_val))


####Doublets identification(单细胞数据可能会有多个细胞在一个液滴，为了减少这个误差进行的操作)####
#要注意：如果做分化发育，很有可能会去掉一些应该存在的细胞，因此可以先保存一个未进行doublets的数据
library(DoubletFinder)
SRR780_val@meta.data[1:5,]
#saveRDS(SRR780_val,file="./SRR780.rds")
#load("./SRR780.rds")
#SRR780_val <- SRR780
SRR780_val_db<-SRR780_val

#先进行降维分析，随机挑选细胞生成假的双胞，计算假双胞与真的分出来的单细胞之间的关系
sweep.res.list_SRR780_val <- paramSweep_v3(SRR780_val_db, PCs = 1:30)
sweep.stats_SRR780_val <- summarizeSweep(sweep.res.list_SRR780_val, GT = FALSE)
bcmvn_SRR780_val<-find.pK(sweep.stats_SRR780_val)
pK_value <- as.numeric(as.character(bcmvn_SRR780_val$pK[bcmvn_SRR780_val$BCmetric == max(bcmvn_SRR780_val$BCmetric)]))
annotations <- SRR780_val_db@meta.data$RNA_snn_res.0.5
homotypic.prop <- modelHomotypic(annotations)  ####get the estimated number
nExp_poi <- round(homotypic.prop*length(SRR780_val_db@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
options(future.globals.maxSize= 50000*1024^2)
#memory.limit(2000)
SRR780_val_db <- doubletFinder_v3(SRR780_val_db, PCs = 1:30, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
SRR780_val_db <- doubletFinder(SRR780_val_db, PCs = 1:30, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
SRR780_val_db@meta.data[1:5,]
table(SRR780_val_db$DF.classifications_0.25_0.01_3563)

SRR780_val$doublets<-SRR780_val_db$DF.classifications_0.25_0.01_3563
SRR780_val@meta.data[1:5,]
####saveRDS###
saveRDS(SRR780_val,file = "./SRR780_val.rds")
DimPlot(SRR780_val,reduction = "tsne",group.by = "doublets")
DimPlot(SRR780_val,reduction = "tsne",group.by = "cluster",label = T)
DimPlot(SRR780_val,reduction = "tsne",label = T)
save.image(file="./tumor.RData")
