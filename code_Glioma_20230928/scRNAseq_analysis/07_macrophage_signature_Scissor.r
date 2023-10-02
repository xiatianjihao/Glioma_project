library(Scissor)
library(ggplot2)

#MDM
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/Scissor/MDM")
sc <- readRDS("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/RDS_file/MDM.rds")
counts <- sc@assays$RNA@counts

sc_dataset <- Seurat_preprocessing(counts, verbose = F)

pdf("sc_umap.pdf")
DimPlot(sc, reduction = 'umap', label = T, label.size = 5)
dev.off()

dir <- "/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database_for_Scissor"

file <- "CGGA_mRNAseq_db_result.xls"
name <- gsub(".xls","",file)
db <- read.csv(paste0(dir,"/",file),header=T,sep="\t",row.names=1)
length(colnames(db))
bulk_dataset <- t(db[,c(3:length(colnames(db)))])
bulk_survival <- db[,c(1,2)]
phenotype <- bulk_survival
colnames(phenotype) <- c("time", "status")
head(phenotype)
infos3 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = 0.01, cutoff = 0.2, family = "cox", Save_file = paste0('Scissor_',name,'_survival.RData'))
Scissor_select <- rep("Background", ncol(sc))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos3$Scissor_pos] <- "Scissor+"
Scissor_select[infos3$Scissor_neg] <- "Scissor-"
sc@meta.data$Scissor <- "NA"
sc@meta.data[names(Scissor_select),]$Scissor <- Scissor_select
p2 <- DimPlot(sc, reduction = 'umap', label = T, label.size = 5,group.by = 'Scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c("Scissor-","Scissor+"))+ggtitle(name)
save(infos3,file=paste0(name,"_infos3.rda"))

name="CGGA_mRNAseq"
Scissor_select <- rep("Background", ncol(sc))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos3$Scissor_pos] <- "Worse survival"
Scissor_select[infos3$Scissor_neg] <- "Better survival"
sc@meta.data$Scissor <- "NA"
sc@meta.data[names(Scissor_select),]$Scissor <- Scissor_select
p2 <- DimPlot(sc, reduction = 'umap', label = F, group.by = 'Scissor', cols = c('grey','indianred1','royalblue'), pt.size = 0.1, order = c("Worse survival","Background"))+xlab("UMAP-1")+ylab("UMAP-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file=paste0(name,"_sc_umap_infos3.pdf"),width=4.5,height=3)

load("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/Scissor/MDM/CGGA_mRNAseq_db_result_infos3.rda")
load("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/Scissor/MDM/Scissor_CGGA_mRNAseq_db_result_survival.RData")
numbers <- length(infos3$Scissor_pos)+length(infos3$Scissor_neg)
result <- reliability.test(X,Y,network,alpha=0.01,family="cox",cell_num=numbers,n=30,nfold=10)
save(result,file="MDM.reliability.rda")
#[1] "Test statistic = 0.570"
#[1] "Reliability significance test p = 0.000"
evaluate_summary <- evaluate.cell("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/Scissor/MDM/Scissor_CGGA_mRNAseq_db_result_survival.RData",infos3,FDR=0.05,bootstrap_n=100)
save(evaluate_summary,file='MDM.evaluate_summary.rda')
evaluate_summary[1:5,1:4]
all(evaluate_summary$`Mean correlation` & as.numeric(gsub('%','',evaluate_summary$`Correlation > 0`)) > 50)
#[1] TRUE

#Microglia
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/Scissor/Microglia")
sc <- readRDS("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/RDS_file/MG.rds")
counts <- sc@assays$RNA@counts

sc_dataset <- Seurat_preprocessing(counts, verbose = F)

pdf("sc_umap.pdf")
DimPlot(sc, reduction = 'umap', label = T, label.size = 5)
dev.off()

dir <- "/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database_for_Scissor"

file <- "CGGA_mRNAseq_db_result.xls"
name <- gsub(".xls","",file)
db <- read.csv(paste0(dir,"/",file),header=T,sep="\t",row.names=1)
length(colnames(db))
bulk_dataset <- t(db[,c(3:length(colnames(db)))])
bulk_survival <- db[,c(1,2)]
phenotype <- bulk_survival
colnames(phenotype) <- c("time", "status")
head(phenotype)
infos3 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = 0.01, cutoff = 0.2, family = "cox", Save_file = paste0('Scissor_',name,'_survival.RData'))
Scissor_select <- rep("Background", ncol(sc))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos3$Scissor_pos] <- "Scissor+"
Scissor_select[infos3$Scissor_neg] <- "Scissor-"
sc@meta.data$Scissor <- "NA"
sc@meta.data[names(Scissor_select),]$Scissor <- Scissor_select
p2 <- DimPlot(sc, reduction = 'umap', label = T, label.size = 5,group.by = 'Scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c("Scissor-","Scissor+"))+ggtitle(name)
save(infos3,file=paste0(name,"_infos3.rda"))

name="CGGA_mRNAseq"
Scissor_select <- rep("Background", ncol(sc))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos3$Scissor_pos] <- "Worse survival"
Scissor_select[infos3$Scissor_neg] <- "Better survival"
sc@meta.data$Scissor <- "NA"
sc@meta.data[names(Scissor_select),]$Scissor <- Scissor_select
p2 <- DimPlot(sc, reduction = 'umap', label = F, group.by = 'Scissor', cols = c('royalblue','grey','indianred1'), pt.size = 0.1, order = c("Worse survival","Background"))+xlab("UMAP-1")+ylab("UMAP-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file=paste0(name,"_sc_umap_infos3.pdf"),width=4.5,height=3)

load("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/Scissor/Microglia/CGGA_mRNAseq_db_result_infos3.rda")
load("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/Scissor/Microglia/Scissor_CGGA_mRNAseq_db_result_survival.RData")
numbers <- length(infos3$Scissor_pos)+length(infos3$Scissor_neg)
result <- reliability.test(X,Y,network,alpha=0.01,family="cox",cell_num=numbers,n=30,nfold=10)
save(result,file="Microglia.reliability.rda")
#[1] "Test statistic = 0.566"
#[1] "Reliability significance test p = 0.033"
evaluate_summary <- evaluate.cell("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/Scissor/Microglia/Scissor_CGGA_mRNAseq_db_result_survival.RData",infos3,FDR=0.05,bootstrap_n=100)
save(evaluate_summary,file='MG.evaluate_summary.rda')
evaluate_summary[1:5,1:4]
all(evaluate_summary$`Mean correlation` & as.numeric(gsub('%','',evaluate_summary$`Correlation > 0`)) > 50)
#[1] TRUE
