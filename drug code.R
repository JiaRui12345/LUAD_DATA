library(readxl)			
library(impute)			
library(limma)			
library(ggplot2)			
library(ggpubr)			

rt1 <- read_excel("/Volumes/TOSHIBA\ EXT/教学资料/王_肺腺癌单细胞_6/第五节课/drug/DTP_NCI60_ZSCORE.xlsx",skip=7)			
colnames(rt1) <- rt1[1,]			
rt1 <- rt1[-1,-c(67,68)]			
rt1 <- rt1[rt1$`FDA status` %in% c("FDA approved","Clinical trial"),]			
rt1 <- rt1[,-c(1,3:6)]			

rt2 <- read_excel("/Volumes/TOSHIBA\ EXT/教学资料/王_肺腺癌单细胞_6/第五节课/drug/RNA__RNA_seq_composite_expression.xls",skip=9)			
colnames(rt2) <- rt2[1,]			
rt2 <- rt2[-1,-c(2:6)]			
write.table(rt2,'/Volumes/TOSHIBA\ EXT/教学资料/王_肺腺癌单细胞_6/第五节课/drug/geneExp.txt',sep='\t',row.names=F,quote=F)			

rt <- as.matrix(rt1)			
rownames(rt) <- rt[,1]			
drug <- rt[,2:ncol(rt)]			
dimnames <- list(rownames(drug),colnames(drug))			
data <- matrix(as.numeric(as.matrix(drug)),nrow=nrow(drug),dimnames=dimnames)			

data <- data[, which(colMeans(!is.na(data)) > 0.5)]			
mat <- impute.knn(data)			
drug <- mat$data			
drug <- avereps(drug)			
write.table(drug,"/Volumes/TOSHIBA\ EXT/教学资料/王_肺腺癌单细胞_6/第五节课/drug/drug.txt",quote=F,sep='\t')			

exp <- read.table("/Volumes/TOSHIBA\ EXT/教学资料/王_肺腺癌单细胞_6/第五节课/drug/geneExp.txt",sep='\t',row.names=1,header=T)			
genes <- readRDS("/Volumes/TOSHIBA\ EXT/教学资料/王_肺腺癌单细胞_6/第五节课/modeling/geneCoef.rds")			
genelist <- genes[,1]			
exp <- exp[genelist,]			

drug <- read.table("/Volumes/TOSHIBA\ EXT/教学资料/王_肺腺癌单细胞_6/第五节课/drug/drug.txt",row.names=1,sep='\t',header=T)			
exp <- exp[,colnames(drug)]			
outTab <- data.frame()			
for (Gene in row.names(exp)){			
  x <- as.numeric(exp[Gene,])		
  for (Drug in row.names(drug)){		
    y <- as.numeric(drug[Drug,])	
    corT <- cor.test(x,y,method="pearson")	
    cor <- corT$estimate	
    pvalue <- corT$p.value	
    if(pvalue < 0.01){	
      outVector <- cbind(Gene, Drug,cor, pvalue)
      outTab <- rbind(outTab,outVector)
    }	
  }		
}			

outTab <- outTab[order(as.numeric(as.vector(outTab$pvalue))),]			
write.table(outTab,"/Volumes/TOSHIBA\ EXT/教学资料/王_肺腺癌单细胞_6/第五节课/drug/drugCor.txt",sep='\t',row.names=F,quote=F)			

plotList_1 <- list()			
corPlotNum <- 16			
for (i in 1:corPlotNum){			
  Gene <- outTab[i,1]		
  Drug <- outTab[i,2]		
  x <- as.numeric(exp[Gene,])		
  y <- as.numeric(drug[Drug,])		
  cor <- sprintf("%.03f",as.numeric(outTab[i,3]))		
  pvalue=0		
  if (as.numeric(outTab[i,4])<0.001){		
    pvalue="p<0.001"	
  }else{		
    pvalue=paste0("p=",sprintf("%,03f",as.numeric(outTab[i,4])))	
  }		
  df1 <- as.data.frame(cbind(x,y))		
  p1 <- ggplot(data=df1,aes(x=x,y=y))+		
    geom_point(size=1)+	
    stat_smooth(method="lm",se=FALSE,formula=y~x)+	
    labs(x="Expression",y="IC50",title=paste0(Gene,", ", Drug),subtitle=paste0("Cor=",cor,", p=",pvalue))+	
    theme(axis.sticks=element_blank(),axis.text.y=element_blank(),axis.text.x=element_blank())+	
    theme_bw()	
  plotList_1[[i]]=p1		
}			

nrow <- ceiling(sqrt(corPlotNum))			
ncol <- ceiling(sqrt(corPlotNum))			
ggarrange(plotlist=plotList_1,nrow=nrow,ncol=ncol)			
