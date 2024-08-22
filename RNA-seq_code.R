library(dplyr)
library(rjson)
library(e1071)
library(parallel)
library(preprocessCore)
library(tibble)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(ConsensusClusterPlus)
library(ggfortify)
library(DESeq2)
library(glmnet)
library(timeROC)
library(limma)
library(GEOquery)
setwd("F:/wang/wang")
#source("/Users/macbookpro/Documents/jupyter/wang/input/cibersort.R")
tcga.meta <- fromJSON(file='./metadata.cart.2023-05-21.json')
tcga.data <- as.data.frame(matrix(NA,nrow=length(tcga.meta),ncol=3))
colnames(tcga.data) <- c('sample.name','file.id','file.name')
for (i in c(1:length(tcga.meta))){
  tcga.data[i,'sample.name'] <- tcga.meta[[i]]$associated_entities[[1]]$entity_submitter_id
  tcga.data[i,'file.id'] <- tcga.meta[[i]]$file_id
  tcga.data[i,'file.name'] <- tcga.meta[[i]]$file_name
}
unzip.file <- './gdc_download_20230521_013108.833593'
for (i in c(1:nrow(tcga.data))){
  sample.count <- read.table(paste(unzip.file,tcga.data[i,'file.id'],tcga.data[i,'file.name'],sep='/'),skip=6,header=F,row.names=1)
  if (i == 1){
    tcga.count <- as.data.frame(sample.count$V8)
  }else{
    tcga.count <- cbind(tcga.count,as.data.frame(sample.count$V8))
  }
}
rownames(tcga.count) <- rownames(sample.count)
colnames(tcga.count) <- tcga.data[,1]
gene.anno <- sample.count[,c(1,2)]
tcga.count <- cbind(gene.anno,tcga.count)
fpkm.id.data <- tcga.count %>% group_by(V2) %>% summarise_at(vars(matches("TCGA")), mean) %>% as.data.frame()
rownames(fpkm.id.data) <- fpkm.id.data$V2
fpkm.id.data$V2 <- NULL
tumor.sample <- c()
tumor.table <- data.frame()
for (i in colnames(fpkm.id.data)){
  tmp.list <- strsplit(i,'-')
  if (startsWith(tmp.list[[1]][4],'0')){
    tumor.sample <- c(tumor.sample,i)
    tumor.table <- rbind(tumor.table,data.frame(ori=i,name=paste(tmp.list[[1]][1],tmp.list[[1]][2],tmp.list[[1]][3],sep='-')))
  }
}
fpkm.id.data <- fpkm.id.data[,tumor.sample]

tmp <- tumor.table %>%
  group_by(name) %>%
  filter(n() == 1) %>%
  ungroup()
fpkm.id.data <- fpkm.id.data[,tmp$ori]
colnames(fpkm.id.data) <- tmp$name
saveRDS(fpkm.id.data,"./fpkm_count.rds")

setwd("F:/LUAD-CAF/step4_cibersortx")
marker <- readRDS("./fb_marker.rds")
union_gene <- intersect(rownames(marker),rownames(fpkm.id.data))
marker <- marker[union_gene,]
CoreAlg <- function(X, y, absolute, abs_method){
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  if(!absolute || abs_method == 'sig.score') w <- (q/sum(q)) #relative space (returns fractions)
  if(absolute && abs_method == 'no.sumto1') w <- q #absolute space (returns scores)
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}

doPerm <- function(perm, X, Y, absolute, abs_method){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while(itor <= perm){
    #print(itor)
    
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr, absolute, abs_method)
    
    mix_r <- result$mix_r
    
    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}
    
    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

CIBERSORT <- function(sig_matrix = lm22, mixture_file, perm, QN = TRUE, absolute, abs_method='sig.score'){
  
  
  if (length(intersect(rownames(mixture_file), rownames(sig_matrix))) == 0){
    stop("None identical gene between eset and reference had been found.
         Check your eset using: intersect(rownames(eset), rownames(reference))")
  }
  
  if(absolute && abs_method != 'no.sumto1' && abs_method != 'sig.score')
    stop("abs_method must be set to either 'sig.score' or 'no.sumto1'")
  
  #read in data
  X<- sig_matrix
  # X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
  Y <- rownames_to_column(mixture_file,var = "symbol")
  #to prevent crashing on duplicated gene symbols, add unique numbers to identical names
  dups <- dim(Y)[1] - length(unique(Y[,1]))
  if(dups > 0) {
    warning(paste(dups," duplicated gene symbol(s) found in mixture file!",sep=""))
    rownames(Y) <- make.names(Y[,1], unique=TRUE)
  }else {rownames(Y) <- Y[,1]}
  Y <- Y[,-1]
  ###################################
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]
  
  P <- perm #number of permutations
  
  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}
  
  #quantile normalization of mixture file
  # library(preprocessCore)
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  #store original mixtures
  Yorig <- Y
  Ymedian <- max(median(Yorig),1)
  
  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  
  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))
  
  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y, absolute, abs_method)$dist)}
  
  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  if(absolute) header <- c(header, paste('Absolute score (',abs_method,')',sep=""))
  
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
  #iterate through mixtures
  while(itor <= mixtures){
    
    y <- Y[,itor]
    
    #standardize mixture
    y <- (y - mean(y)) / sd(y)
    
    #run SVR core algorithm
    result <- CoreAlg(X, y, absolute, abs_method)
    
    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    
    if(absolute && abs_method == 'sig.score') {
      w <- w * median(Y[,itor]) / Ymedian
    }
    
    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
    
    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(absolute) out <- c(out, sum(w))
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}
    
    itor <- itor + 1
    
  }
  
  #save results
  # write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
  
  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  if(!absolute){colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")}
  else{colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE",paste('Absolute score (',abs_method,')',sep=""))}
  obj
}
result <- CIBERSORT(marker,fpkm.id.data,perm=1000,QN=T,absolute=FALSE)
saveRDS(result,"./cibersort.rds")


setwd("F:/wang/wang/input")
tcga.clinical <- read.csv('./clinical.tsv',header=T,sep='\t')
tcga.clinical <- dplyr::select(tcga.clinical,c('case_submitter_id','gender','age_at_index','vital_status','days_to_death','days_to_last_follow_up'))
tcga.clinical <- tcga.clinical[duplicated(tcga.clinical),]
rownames(tcga.clinical) <- tcga.clinical$case_submitter_id
for (i in rownames(tcga.clinical)){
  if (tcga.clinical[i,"days_to_death"]=="'--"){
    tcga.clinical[i,"days_to_death"] <- tcga.clinical[i,"days_to_last_follow_up"]
  }
}
tcga.clinical$days_to_last_follow_up <- NULL
tcga.clinical <- filter(tcga.clinical,days_to_death != "'--")
tcga.clinical$days_to_death <- as.numeric(tcga.clinical$days_to_death)/365
tcga.clinical$vital_status <- ifelse(tcga.clinical$vital_status=='Alive',0,1)
union.sample <- intersect(colnames(fpkm.id.data),rownames(tcga.clinical))
fpkm.id.data <- fpkm.id.data[,union.sample]
tcga.clinical <- tcga.clinical[union.sample,]
#去除小于30天的
tcga.clinical <- filter(tcga.clinical,days_to_death>30/365)
fpkm.id.data <- fpkm.id.data[,rownames(tcga.clinical)]

union.sample <- intersect(colnames(fpkm.id.data),rownames(tcga.clinical))
fpkm.id.data <- fpkm.id.data[,union.sample]
tcga.clinical <- tcga.clinical[union.sample,]


result <- as.data.frame(result[union.sample,c(1:7)])
tcga.clinical$immune <- result$Immune_related.CAF/rowSums(result)
tcga.clinical$MCAF <- result$MyCAF/rowSums(result)
ggdensity(tcga.clinical, x = "MCAF", 
          fill = "#0073C2FF", color = "#0073C2FF",
          add = "median", rug = TRUE)
tcga.clinical$cell_group <- ifelse(tcga.clinical$MCAF>median(tcga.clinical$MCAF),"high",'low')
diff <- survdiff(Surv(days_to_death,as.numeric(vital_status))~cell_group,data=tcga.clinical)
p.value <- 1-pchisq(diff$chisq,df=0.5)
p.value <- 1-pchisq(diff$chisq,df=length(diff$n)-1)
p.value <- signif(p.value,4)
p.value <- format(p.value,scientific=T)
fit <- survfit(Surv(days_to_death,as.numeric(vital_status))~cell_group, data=tcga.clinical)
ggsurvplot(fit,data=tcga.clinical,conf.int=T,pval=paste0("p=",p.value),
           pval.size=5,legend.title="Risk",
           legend.labs=c("High MCAF", "Low MCAF"),
           xlab="Time(years)",break.time.by=1,
           palette=c("red","blue"),
           risk.table=F,risk.table.title="",risk.table.height=.25)
ggdensity(tcga.clinical, x = "immune", 
          fill = "#0073C2FF", color = "#0073C2FF",
          add = "median", rug = TRUE)
tcga.clinical$cell_group_immune <- ifelse(tcga.clinical$immune>median(tcga.clinical$immune),"high",'low')
diff <- survdiff(Surv(days_to_death,as.numeric(vital_status))~cell_group_immune,data=tcga.clinical)
p.value <- 1-pchisq(diff$chisq,df=0.5)
p.value <- 1-pchisq(diff$chisq,df=length(diff$n)-1)
p.value <- signif(p.value,4)
p.value <- format(p.value,scientific=T)
fit <- survfit(Surv(days_to_death,as.numeric(vital_status))~cell_group_immune, data=tcga.clinical)
ggsurvplot(fit,data=tcga.clinical,conf.int=T,pval=paste0("p=",p.value),
           pval.size=5,legend.title="Risk",
           legend.labs=c("High immune", "Low immune"),
           xlab="Time(years)",break.time.by=1,
           palette=c("red","blue"),
           risk.table=F,risk.table.title="",risk.table.height=.25)

#getmarker		
library(Seurat)	
setwd("F:/LUAD-CAF/step4_cibersortx")
seurat.obj <- readRDS("../step1_merge/sce.mergeTEN.rds")		
fibroblast <- readRDS("../step2_caf/fb_seurat.rds")		
fibroblast$celltype <- fibroblast@active.ident		
mcafs <- subset(fibroblast, subset = celltype == "MyCAF")		
mcafs <- colnames(mcafs)		
seurat.obj$celltype <- seurat.obj@active.ident		
tmp <- data.frame(cells=colnames(seurat.obj),type=seurat.obj$celltype)		
tmp$type <- as.character(tmp$type)		
tmp[mcafs,"type"] <- "MyCAF"		
seurat.obj$celltype <- tmp$type		
seurat.obj <- SetIdent(seurat.obj,value='celltype')		
fb.de.markers <- FindMarkers(seurat.obj, ident.1 = "MyCAF")		
mcaf.marker <- fb.de.markers %>% filter(p_val_adj<0.01 & avg_log2FC>2) %>% rownames()		
#mcaf.marker <- intersect(mcaf.marker,rownames(tcga.count))		
mycaf_DGEs <- fb.de.markers[mcaf.marker,]
write.csv(mycaf_DGEs,"./mycaf_DGEs.csv")


#consesus
#setwd("F:/wang/wang")
#mcaf.marker <- readRDS("/Users/macbookpro/Documents/jupyter/wang/input/mcaf_marker.rds")
#mcaf.marker <- intersect(mcaf.marker,rownames(fpkm.id.data))
#cluster.count <- fpkm.id.data[mcaf.marker,]
#cluster.count = sweep(cluster.count,1, apply(cluster.count,1,median,na.rm=T))
#cluster.count<-as.matrix(cluster.count)
#results = ConsensusClusterPlus(cluster.count,maxK=5,reps=50,pItem=0.8,pFeature=1, title='MCAF marker clustering',clusterAlg="hc",distance="pearson",seed=7,plot=NULL)
#tcga.clinical$con <- paste("cluster",results[[2]]$consensusClass,sep='_')
#diff <- survdiff(Surv(days_to_death,as.numeric(vital_status))~con,data=tcga.clinical)
#p.value <- 1-pchisq(diff$chisq,df=0.1)
#p.value <- 1-pchisq(diff$chisq,df=length(diff$n)-1)
#p.value <- signif(p.value,4)
#p.value <- format(p.value,scientific=T)
#fit <- survfit(Surv(days_to_death,as.numeric(vital_status))~con, data=tcga.clinical)
#ggsurvplot(fit,data=tcga.clinical,conf.int=T,pval=paste0("p=",p.value),
#           pval.size=5,legend.title="Risk",
#           legend.labs=c("Cluster1", "Cluster2"),
#           xlab="Time(years)",break.time.by=1,
#           palette=c("red","blue"),
#           risk.table=F,risk.table.title="",risk.table.height=.25)
#pca_res <- prcomp(t(cluster.count), scale. = TRUE)
#autoplot(pca_res,data=tcga.clinical,colour='con',frame=F,type='norm')

#de
unzip.file <- './input/gdc_download_20230521_013108.833593'
for (i in c(1:nrow(tcga.data))){
  sample.count <- read.table(paste(unzip.file,tcga.data[i,'file.id'],tcga.data[i,'file.name'],sep='/'),skip=6,header=F,row.names=1)
  if (i == 1){
    tcga.count <- as.data.frame(sample.count$V4)
  }else{
    tcga.count <- cbind(tcga.count,as.data.frame(sample.count$V4))
  }
}
rownames(tcga.count) <- rownames(sample.count)
colnames(tcga.count) <- tcga.data[,1]
gene.anno <- sample.count[,c(1,2)]
tcga.count <- cbind(gene.anno,tcga.count)
raw.id.data <- tcga.count %>% group_by(V2) %>% summarise_at(vars(matches("TCGA")), max) %>% as.data.frame()
rownames(raw.id.data) <- raw.id.data$V2
raw.id.data$V2 <- NULL

clinical <- data.frame()
for (i in colnames(raw.id.data)){
  tmp.list <- strsplit(i,'-')
  if (startsWith(tmp.list[[1]][4],'0')){
    clinical <- rbind(clinical,data.frame(sample=i,type="tumor"))
  }else{
    clinical <- rbind(clinical,data.frame(sample=i,type="normal"))
  }
}
table(clinical$type)
dds <- DESeqDataSetFromMatrix(countData = raw.id.data,
                              colData = clinical,
                              design = ~ type)
dds <- DESeq(dds)
res <- results(dds)
res <- as.data.frame(res)
res$diffexpressed <- 'Not significant'
genes <- res %>% filter(padj<0.01&log2FoldChange< -1) %>% rownames()
res[genes,"diffexpressed"] <- 'Downregulated'
genes <- res %>% filter(padj<0.01&log2FoldChange> 1) %>% rownames()
res[genes,"diffexpressed"] <- 'Upregulated'
ggplot(data = res, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "black", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) + 
  coord_cartesian(ylim = c(0, 180), xlim = c(-6, 11)) + 
  labs(color = 'Severe', 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"q-value")) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  ggtitle('Cancer tissue VS Normal lung tissue')
diff.gene <- res %>% filter(abs(log2FoldChange)>1&padj<0.01) %>% rownames()
DEGs <- res[diff.gene,]
write.csv(DEGs,"./DEGs.csv")
write.csv(tcga.clinical,"./tcga.clinical.csv")
write.csv(geo.sample,"./geo.sample.csv")

genename <- as.character(rownames(res))
gene_map <- mapIds(org.Hs.eg.db, genename, 'ENTREZID', 'SYMBOL')
geneList <- res$log2FoldChange
names(geneList) <- gene_map
geneList<-na.omit(geneList)
geneList = geneList[order(geneList,decreasing = TRUE)]
#GO gsea
GO <- gseGO(
  geneList, #gene_fc
  ont = "ALL",
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
)
GO_result <- GO@result
write.csv(GO_result,"./GO_result.csv")

gseaplot2(GO, c(1,2,3,29,30), pvalue_table = TRUE)

gseaplot2(GO, c(4:8), pvalue_table = TRUE)

#model
genes <- intersect(mcaf.marker,diff.gene)
tcga.clinical <- tcga.clinical %>% filter(days_to_death>30/365)
x <- as.matrix(t(fpkm.id.data[genes,rownames(tcga.clinical)]))
y <- data.matrix(Surv(as.numeric(tcga.clinical$days_to_death),as.numeric(tcga.clinical$vital_status)))
set.seed(56)
fit <- glmnet(x,y,family='cox',maxit=10000)
set.seed(56)
cvfit <- cv.glmnet(x,y,family="cox",maxit=10000,nfole=5,type.measure='C')
plot(fit)
plot(cvfit)
coef <- coef(fit,s=cvfit$lambda.1se)
index <- which(coef!=0)
actCoef <- coef[index]
lassoGene <- row.names(coef)[index]
geneCoef <- cbind(Gene=lassoGene,Coef=actCoef)
tcga.clinical$risk_score <- 0
for (i in rownames(tcga.clinical)){
  rs = 0
  for (j in c(1:nrow(geneCoef))){
    gene <- geneCoef[j,1]
    rs = rs+as.numeric(geneCoef[j,2])*as.numeric(fpkm.id.data[gene,i])
  }
  tcga.clinical[i,'risk_score'] <- rs
}
tcga.clinical$group <- ifelse(tcga.clinical$risk_score > median(tcga.clinical$risk_score), "high", "low")
rs.curve <- tcga.clinical[order(tcga.clinical$risk_score),]
riskClass <- rs.curve[,'group']
lowLength <- length(riskClass[riskClass=='low'])
highLength <- length(riskClass[riskClass=='high'])
lowMax <- max(rs.curve$risk_score[riskClass=='low'])
line <- rs.curve[,'risk_score']
plot(line,type='p',pch=20,xlab='Patients (increasing risk score)', 
     ylab="Risk score",col=c(rep('blue',lowLength),rep("red",highLength)))
abline(h=lowMax,v=lowLength,lty=2)
legend("topleft",c("High risk","low Risk"),bty='n',pch=19,col=c("red","blue"),cex=1.2)
color <- as.vector(rs.curve$vital_status)
color[color==1] = "red"
color[color==0] = "blue"
plot(as.numeric(rs.curve$days_to_death),pch=19,
     xlab="Patients (increasing risk score)", ylab="Survival time (years)",col=color)
abline(v=lowLength,lty=2)
legend("topleft",c("Dead","Alive"),bty='n',pch=19,col=c("red","blue"),cex=1.2)
diff <- survdiff(Surv(days_to_death,as.numeric(vital_status))~group,data=tcga.clinical)
#p.value <- 1-pchisq(diff$chisq,df=.5)
p.value <- 1-pchisq(diff$chisq,df=length(diff$n)-1)
p.value <- signif(p.value,4)
p.value <- format(p.value,scientific=T)
fit <- survfit(Surv(days_to_death,as.numeric(vital_status))~group, data=tcga.clinical)
ggsurvplot(fit,data=tcga.clinical,conf.int=T,pval=paste0("p=",p.value),
           pval.size=5,legend.title="Risk",
           legend.labs=c("high risk", "low risk"),
           xlab="Time(years)",break.time.by=1,
           palette=c("red","blue"),
           risk.table.title="",risk.table.height=.25,
           risk.table = TRUE)
roc <- timeROC(T=tcga.clinical$days_to_death,delta=tcga.clinical$vital_status,marker=tcga.clinical$risk_score,
               cause=1,weighting='aalen',times=c(.8,2.4,3.3),ROC=T)
plot(roc,time=0.8,col='green',title=F,lwd=2)
plot(roc,time=2.4,col='blue',add=T,title=F,lwd=2)
plot(roc,time=3.3,col='red',add=T,title=F,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",roc$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",roc$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",roc$AUC[3]))),
       col=c("green","blue","red"),lwd=2,bty='n')
clinical <- read.csv('./input/clinical.tsv',header=T,sep='\t')
rownames(clinical) <- make.unique(clinical$case_submitter_id)
clinical <- clinical[rownames(tcga.clinical),]
clinical <- cbind(clinical[,c('age_at_index','gender','ajcc_pathologic_t','ajcc_pathologic_stage')],tcga.clinical[,c('risk_score','days_to_death','vital_status')])
clinical <- subset(clinical,subset=ajcc_pathologic_stage!="'--")
clinical$age_at_index <- as.numeric(clinical$age_at_index)
clinical <- clinical[!is.na(clinical$age_at_index),]
clinical$gender <- ifelse(clinical$gender=='female',0,1)
clinical$ajcc_pathologic_stage <- gsub('Stage IA','1',clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_stage <- gsub('Stage IB','1',clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_stage <- gsub('Stage IIA','2',clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_stage <- gsub('Stage IIB','2',clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_stage <- gsub('Stage IV','4',clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_stage <- gsub('Stage IIIA','3',clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_stage <- gsub('Stage IIIB','3',clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_stage <- gsub('Stage III','3',clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_stage <- gsub('Stage II','2',clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_stage <- gsub('Stage I','1',clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_stage <- as.numeric(clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_t <- gsub('T1a','1',clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_t <- gsub('T1b','1',clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_t <- gsub('T1','1',clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_t <- gsub('T2a','2',clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_t <- gsub('T2b','2',clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_t <- gsub('T2','2',clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_t <- gsub('T3','3',clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_t <- gsub('T4','4',clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_t <- gsub('TX','5',clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_t <- as.numeric(clinical$ajcc_pathologic_t)

uniTab <- data.frame()
for (i in c('age_at_index','gender','ajcc_pathologic_t','ajcc_pathologic_stage','risk_score')){
  cox <- coxph(Surv(days_to_death,vital_status) ~ clinical[,i],data=clinical)
  coxSummary <- summary(cox)
  uniTab <- rbind(uniTab,cbind(
    id=i,
    HR=coxSummary$conf.int[,"exp(coef)"],
    HR.95L=coxSummary$conf.int[,"lower .95"],
    HR.95H=coxSummary$conf.int[,"upper .95"],
    pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}

uniTab$HR.95L <- as.numeric(uniTab$HR.95L)
uniTab$HR.95H <- as.numeric(uniTab$HR.95H)
uniTab$CI <- paste(round(uniTab$HR.95L,2),round(uniTab$HR.95H,2),sep=', ')
uniTab$p.value <- round(as.numeric(uniTab$pvalue),8)
uniTab$index <- c(1:5)
write.csv(uniTab,"./uniTab.csv")
library(forestplot)
BiocManager::install("forestplot")
tabletext<- cbind(c("id",uniTab$id),
                  c("HR",format(round(as.numeric(uniTab$HR),3),nsmall = 3)),
                  c("lower 95%CI",format(round(as.numeric(uniTab$HR.95L),3),nsmall = 3)),
                  c("upper 95%CI",format(round(as.numeric(uniTab$HR.95H),3),nsmall = 3)),
                  c("pvalue",format(round(as.numeric(uniTab$p.value),3),nsmall = 3)))
##3.2 绘制森林图
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(uniTab$HR)),
           lower=c(NA,as.numeric(uniTab$HR.95L)), 
           upper=c(NA,as.numeric(uniTab$HR.95H)),
           graph.pos=5,# 图在表中的列位置
           graphwidth = unit(.25,"npc"),# 图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",# box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),# box颜色
           
           boxsize=0.4,# box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=T,# 显示区间
           zero=1,# zero线横坐标
           lwd.zero=1.5,# zero线宽
           xticks = c(1,1.5,2,2.5,3,3.5,4,4.5,5),# 横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1.2),# 各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"), # 在第一行上面画黑色实线
                           "2" = gpar(lwd=1.5, col="black"), # 在第一行标题行下画黑色实线
                           "7" = gpar(lwd=2, col="black")), # 在最后一行上画黑色实线
           lineheight = unit(.75,"cm"),# 固定行高
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
dev.off()


#vil
y.raw <- read.delim("./input/GSE42127_non-normalized.txt.gz",row.names=NULL)
row.names(y.raw) <- make.unique(y.raw$ID_REF)
y.raw$ID_REF <- NULL
y.bgcorrected <- backgroundCorrect(y.raw, method="normexp")
y.norm <- normalizeBetweenArrays(log2(y.bgcorrected), method="quantile")
anno <- read.table("./input/GPL6884-11607.txt", header=T, comment.char="#", sep='\t', fill=NA)
anno <- anno %>% filter(ILMN_Gene %in% geneCoef[,1])

predict <- y.norm[anno$ID,]
predict <- cbind(data.frame(gene=anno[,"ILMN_Gene"]),predict)
predict <- predict %>% group_by(gene) %>% summarise_at(vars(colnames(predict)[2:177]), mean) %>% as.data.frame()
rownames(predict) <- predict$gene
predict$gene <- NULL
eSet <- getGEO(file='./input/GSE42127_series_matrix.txt.gz',AnnotGPL = F,getGPL = F)
geo.sample <- pData(eSet)
geo.sample <- subset(geo.sample,subset=`histology:ch1`=='Adenocarcinoma')
#geo.sample <- geo.sample %>% select(c("description",'overall survival months:ch1','survival status:ch1','histology:ch1'))
geo.sample <- dplyr::select(geo.sample,c("description",'overall survival months:ch1','survival status:ch1','histology:ch1'))
colnames(geo.sample) <- c("sample",'days_to_death','vital_status','type')
geo.sample$days_to_death <- as.numeric(geo.sample$days_to_death)/12
geo.sample$vital_status <- ifelse(geo.sample$vital_status=='A',0,1)
geo.sample$sample <- paste0('X',gsub('-','.',geo.sample$sample))
rownames(geo.sample) <- geo.sample$sample
geo.sample$risk_score <- 0
rownames(geneCoef) <- geneCoef[,1]
for (i in rownames(geo.sample)){
  rs = 0
  for (j in c(1:nrow(predict))){
    gene <- rownames(predict)[j]
    rs = rs+as.numeric(geneCoef[gene,2])*as.numeric(predict[gene,i])
  }
  geo.sample[i,'risk_score'] <- rs
}
geo.sample$group <- ifelse(geo.sample$risk_score > median(geo.sample$risk_score), "high", "low")
rs.curve <- geo.sample[order(geo.sample$risk_score),]
riskClass <- rs.curve[,'group']
lowLength <- length(riskClass[riskClass=='low'])
highLength <- length(riskClass[riskClass=='high'])
lowMax <- max(rs.curve$risk_score[riskClass=='low'])
line <- rs.curve[,'risk_score']
plot(line,type='p',pch=20,xlab='Patients (increasing risk score)', 
     ylab="Risk score",col=c(rep('blue',lowLength),rep("red",highLength)))
abline(h=lowMax,v=lowLength,lty=2)
legend("topleft",c("High risk","low Risk"),bty='n',pch=19,col=c("red","blue"),cex=1.2)
color <- as.vector(rs.curve$vital_status)
color[color==1] = "red"
color[color==0] = "blue"
plot(as.numeric(rs.curve$days_to_death),pch=19,
     xlab="Patients (increasing risk score)", ylab="Survival time (years)",col=color)
abline(v=lowLength,lty=2)
legend("topleft",c("Dead","Alive"),bty='n',pch=19,col=c("red","blue"),cex=1.2)
diff <- survdiff(Surv(days_to_death,as.numeric(vital_status))~group,data=geo.sample)
#p.value <- 1-pchisq(diff$chisq,df=.5)
p.value <- 1-pchisq(diff$chisq,df=length(diff$n)-1)
p.value <- signif(p.value,4)
p.value <- format(p.value,scientific=T)
fit <- survfit(Surv(days_to_death,as.numeric(vital_status))~group, data=geo.sample)
ggsurvplot(fit,data=geo.sample,conf.int=T,pval=paste0("p=",p.value),
           pval.size=5,legend.title="Risk",
           legend.labs=c("high risk", "low risk"),
           xlab="Time(years)",break.time.by=1,
           palette=c("red","blue"),
           risk.table.title="",risk.table.height=.25,
           risk.table = TRUE)
roc <- timeROC(T=geo.sample$days_to_death,delta=geo.sample$vital_status,marker=geo.sample$risk_score,
               cause=1,weighting='aalen',times=c(1.2,2.2,2.8),ROC=T)
plot(roc,time=1.2,col='green',title=F,lwd=2)
plot(roc,time=2.2,col='blue',add=T,title=F,lwd=2)
plot(roc,time=2.8,col='red',add=T,title=F,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",roc$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",roc$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",roc$AUC[3]))),
       col=c("green","blue","red"),lwd=2,bty='n')
