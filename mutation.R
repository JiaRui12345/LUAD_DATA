#load package
library(rjson)
library(maftools)
library(ggpubr)
BiocManager::install("maftools")
#load file
setwd("F:/LUAD-CAF")
tcga.meta <- fromJSON(file='./RAW/TCGA_MUTA/metadata.cart.2023-06-28.json')
"tcga.data <- as.data.frame(matrix(NA,nrow=length(tcga.meta),ncol=4))"
"colnames(tcga.data) <- c('sample.name','sample','file.id','file.name')"
for (i in c(1:length(tcga.meta))){
"  tcga.data[i,'sample.name'] <- tcga.meta[[i]]$associated_entities[[1]]$entity_submitter_id"
"  tcga.data[i,'file.id'] <- tcga.meta[[i]]$file_id"
"  tcga.data[i,'file.name'] <- tcga.meta[[i]]$file_name"
"  tcga.data[i,'sample'] <- paste(unlist(strsplit(tcga.meta[[i]]$associated_entities[[1]]$entity_submitter_id,split='-')[[1]][1:3]),collapse='-')"
}
#survival.data <- readRDS("/tcga_group.rds")
"survival.data <- read.csv(""./step5_bulk consensus/wang/wang/tcga.clinical1.csv"",row.names=1)"
geneCoef <- read.csv("./step5_bulk consensus/wang/wang/lasso_coef.txt.csv")
"genes <- geneCoef[,1]"

#low risk
low.sample <- survival.data %>% filter(group == "low") %>% rownames()
"low.file <- filter(tcga.data,sample %in% low.sample)"
rownames(low.file) <- make.unique(low.file$sample)
"low.file <- low.file[intersect(tcga.data$sample,low.sample),]"
all_mut <- data.frame()
for (i in 1:nrow(low.file)) {
"  file <- paste(""./RAW/TCGA_MUTA/gdc_download_20230628_114100.554550"","
"                low.file[i,""file.id""],low.file[i,""file.name""],sep='/')"
"  mut <- read.delim(file,skip = 7, header = T, fill = TRUE,sep = ""\t"")"
"  all_mut <- rbind(all_mut,mut)"
}
all_mut <- read.maf(all_mut)
tmb_table1 = tmb(maf = all_mut)
"write.csv(tmb_table1,""./step6_mutation_ppi_drug/tmb_table_low.csv"")"
dev.off()
"oncoplot(maf = all_mut,"
"         top = 10, "
"         fontSize = 0.6, "
         showTumorSampleBarcodes = F)
"oncoplot(maf = all_mut,"
"         genes=genes, "
"         fontSize = 0.6, "
         showTumorSampleBarcodes = F)

#high risk
high.sample <- survival.data %>% filter(group=='high') %>% rownames()
"high.file <- filter(tcga.data,sample %in% high.sample)"
rownames(high.file) <- make.unique(high.file$sample)
"high.file <- high.file[intersect(tcga.data$sample,high.sample),]"
all_mut <- data.frame()
for (i in 1:nrow(high.file)) {
"  file <- paste(""./RAW/TCGA_MUTA/gdc_download_20230628_114100.554550"","
"                high.file[i,""file.id""],high.file[i,""file.name""],sep='/')"
"  mut <- read.delim(file,skip = 7, header = T, fill = TRUE,sep = ""\t"")"
"  all_mut <- rbind(all_mut,mut)"
}
all_mut <- read.maf(all_mut)
dev.off()
tmb_table2 = tmb(maf = all_mut)
"write.csv(tmb_table2,""./step6_mutation_ppi_drug/tmb_table_high.csv"")"
#tmb_table <- as.data.frame(tmb_table)
"oncoplot(maf = all_mut,"
"         top = 10, "
"         fontSize = 0.6, "
         showTumorSampleBarcodes = F)
"oncoplot(maf = all_mut,"
"         genes=genes, "
"         fontSize = 0.6, "
         showTumorSampleBarcodes = F)

#tmp
#tmb
"all.sample <- intersect(rownames(survival.data),tcga.data$sample)"
rownames(tcga.data) <- make.unique(tcga.data$sample)
"tcga.data <- tcga.data[all.sample,]"
all_mut <- data.frame()
for (i in 1:nrow(tcga.data)) {
"  file <- paste(""./RAW/TCGA_MUTA/gdc_download_20230628_114100.554550"","
"                tcga.data[i,""file.id""],tcga.data[i,""file.name""],sep='/')"
"  mut <- read.delim(file,skip = 7, header = T, fill = TRUE,sep = ""\t"")"
"  all_mut <- rbind(all_mut,mut)"
}
all_mut <- read.maf(all_mut)
tmb_table = tmb(maf = all_mut)
tmb_table <- as.data.frame(tmb_table)
dev.off()
"write.csv(tmb_table_all,""./step6_mutation_ppi_drug/tmb_table_all.csv"")"
tmb_table$sample <- "NA"
for (i in 1:nrow(tmb_table)){
"  tmb_table[i,""sample""] <- paste(unlist(strsplit(as.character(tmb_table[i,1]),split='-')[[1]][1:3]),collapse='-')"
}
"tmb.plot <- cbind(tmb_table[,c(""total_perMB"",""sample"")],survival.data[tmb_table$sample,c(""risk_score"",""group"")])"
"sp <- ggscatter(tmb.plot, x = ""risk_score"", y = ""total_perMB"","
"                add = ""reg.line"",  # Add regressin line"
"                add.params = list(color = ""blue"", fill = ""lightgray""), # Customize reg. line"
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
"sp + stat_cor(method = ""pearson"", label.x = 1.5, label.y = 20)"
dev.off()
