library(GSVA)
library(GSEABase)
library(RColorBrewer)
library(pheatmap)
library(survival)
library(survminer)
library(grid)
library(gridExtra)


setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/2023_revision/single_gene")

##### TCGA GBMLGG
file <- "/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database/download_GlioVis/2023-07-12_TCGA_GBMLGG_expression.txt"
pheno_file <- "/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database/download_GlioVis/2023-07-12_TCGA_GBMLGG_pheno.txt"
pheno <- read.csv(pheno_file,header=T,sep="\t",row.names=1)
#pheno <- pheno[which(pheno$Histology == "GBM"),]
DB <- read.csv(file,header=T,sep="\t",row.names=1)
commP <- intersect(rownames(pheno),rownames(DB))
DB <- DB[commP,]
pheno <- pheno[commP,]
colnames(DB) <- gsub("HLA\\.","HLA-",colnames(DB))
pheno[which(is.na(pheno$Grade)),"Grade"] <- "no"
pheno[which(is.na(pheno$IDH.status)),"IDH.status"] <- "no"
Grade <- pheno[rownames(DB),"Grade"]
IDH.status <- pheno[rownames(DB),"IDH.status"]
Time <- pheno[rownames(DB),"survival"]
Exist <- pheno[rownames(DB),"status"]
DB2 <- data.frame(Time=Time,Exist=Exist,Grade=Grade,IDH.status=IDH.status,DB[,c(1:length(colnames(DB)))])
DB3 <- DB2[which(DB2$IDH.status == "WT"),]
DB <- DB3[,-c(3,4)]

gene <- "AHR"
dbname <- "TCGA_GBM_LGG_IDHWT"
xsit=75
info <- DB[,c("Time","Exist",gene)]
deg.data <- info
deg.data$group <- ""

deg.data$group[which(deg.data[,gene] >= median(info[,gene]))] <- paste0(gene,"-high")
deg.data$group[which(deg.data[,gene] < median(info[,gene]))] <- paste0(gene,"-low")

dat <- deg.data[which(deg.data$group != ""),]
group <- dat$group
sfit <- survfit(Surv(Time, Exist) ~ group, data = dat)
pdf(paste0(dbname,"-",gene,"-survival.pdf"),width=4,height=5)
p <- ggsurvplot(sfit, pval = TRUE, risk.table = TRUE,risk.table.y.text = FALSE,risk.table.col = "strata", linetype = "solid", pval.coord = c(xsit, 1), font.x = c(12, face = "bold"),font.y = c(12, face = "bold"),  xlab = "Time (Months)", ggtheme = theme_survminer(), legend.labs=c(paste0(gene,"-high"),paste0(gene,"-low")),legend.title="", palette = c("#EE0000FF","#3B4992FF"))+guides(colour = guide_legend(nrow = 2),risk.table.fontsize=8)+guides(colour = guide_legend(nrow = 2))
p$table <- p$table + theme(axis.line = element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
print(p)
dev.off()


#CGGA IDH-WT
file="/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database_CGGA/database/CGGA_mRNAseq_db_all_glioma.txt"

mRNAseq_325_info <- read.csv("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database_CGGA/ori_data/CGGA.mRNAseq_325_clinical.20200506.xls",header=T,sep="\t")
mRNAseq_693_info <- read.csv("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database_CGGA/ori_data/CGGA.mRNAseq_693_clinical.20200506.xls",header=T,sep="\t")
pheno <- rbind(mRNAseq_325_info,mRNAseq_693_info)
rownames(pheno) <- pheno$CGGA_ID

DB <- read.csv(file,header=T,sep="\t",row.names=1)[rownames(pheno),]
commP <- intersect(rownames(pheno),rownames(DB))
DB <- DB[commP,]
pheno <- pheno[commP,]
colnames(DB) <- gsub("HLA\\.","HLA-",colnames(DB))
pheno[which(is.na(pheno$Grade)),"Grade"] <- "no"
pheno[which(is.na(pheno$IDH_mutation_status)),"IDH_mutation_status"] <- "no"
Grade <- pheno[rownames(DB),"Grade"]
IDH.status <- pheno[rownames(DB),"IDH_mutation_status"]
DB2 <- data.frame(DB[,c(1:2)],Grade=Grade,IDH.status=IDH.status,DB[,c(3:length(colnames(DB)))])
DB2[,"Time"] <- DB2[,"Time"]/30
DB3 <- DB2[which(DB2$IDH.status == "Wildtype"),]
DB <- DB3[,-c(3,4)]

gene <- "AHR"
dbname <- "CGGA_allglioma_IDHWT"
xsit=50
info <- DB[,c("Time","Exist",gene)]
deg.data <- info
deg.data$group <- ""

deg.data$group[which(deg.data[,gene] >= median(info[,gene]))] <- paste0(gene,"-high")
deg.data$group[which(deg.data[,gene] < median(info[,gene]))] <- paste0(gene,"-low")

dat <- deg.data[which(deg.data$group != ""),]
group <- dat$group
sfit <- survfit(Surv(Time, Exist) ~ group, data = dat)
pdf(paste0(dbname,"-",gene,"-survival.pdf"),width=4,height=5)
p <- ggsurvplot(sfit, pval = TRUE, risk.table = TRUE,risk.table.y.text = FALSE,risk.table.col = "strata", linetype = "solid", pval.coord = c(xsit, 1), font.x = c(12, face = "bold"),font.y = c(12, face = "bold"),  xlab = "Time (Months)", ggtheme = theme_survminer(), legend.labs=c(paste0(gene,"-high"),paste0(gene,"-low")),legend.title="", palette = c("#EE0000FF","#3B4992FF"))+guides(colour = guide_legend(nrow = 2),risk.table.fontsize=8)+guides(colour = guide_legend(nrow = 2))
p$table <- p$table + theme(axis.line = element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
print(p)
dev.off()


##########################two genes##########################
#two gene survival
# TCGA IDHWT
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/2023_revision/two_genes")
genes1 <- c("SPP1","SPP1","B2M","CCL3","MIF")
genes2 <- c("CD44","ITGA4","KLRC1","CCR5","CD74")

for(i in c(1:5)){
	gene1 <- genes1[i]
	gene2 <- genes2[i]
	file <- "/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database/download_GlioVis/2023-07-12_TCGA_GBMLGG_expression.txt"
	pheno_file <- "/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database/download_GlioVis/2023-07-12_TCGA_GBMLGG_pheno.txt"
	pheno <- read.csv(pheno_file,header=T,sep="\t",row.names=1)
	#pheno <- pheno[which(pheno$Histology == "GBM"),]
	DB <- read.csv(file,header=T,sep="\t",row.names=1)
	commP <- intersect(rownames(pheno),rownames(DB))
	DB <- DB[commP,]
	pheno <- pheno[commP,]
	colnames(DB) <- gsub("HLA\\.","HLA-",colnames(DB))
	pheno[which(is.na(pheno$Grade)),"Grade"] <- "no"
	pheno[which(is.na(pheno$IDH.status)),"IDH.status"] <- "no"
	Grade <- pheno[rownames(DB),"Grade"]
	IDH.status <- pheno[rownames(DB),"IDH.status"]
	Time <- pheno[rownames(DB),"survival"]
	Exist <- pheno[rownames(DB),"status"]
	DB2 <- data.frame(Time=Time,Exist=Exist,Grade=Grade,IDH.status=IDH.status,DB[,c(1:length(colnames(DB)))])
	DB3 <- DB2[which(DB2$IDH.status == "WT"),]
	DB <- DB3[,-c(3,4)]
	colnames(DB) <- gsub("HLA\\.","HLA-",colnames(DB))
	dbname <- "TCGA_GBMLGG_IDHWT"
	xsit=50
	info <- DB[,c("Time","Exist",gene1,gene2)]
	deg.data <- info
	deg.data$group <- ""
	deg.data$group[which(deg.data[,gene1] >= median(info[,gene1]) & (deg.data[,gene2] >= median(info[,gene2])))] <- paste0(gene1,"-high + ",gene2,"-high")
	deg.data$group[which(deg.data[,gene1] < median(info[,gene1]) & (deg.data[,gene2] < median(info[,gene2])))] <- paste0(gene1,"-low + ",gene2,"-low")
	dat <- deg.data[which(deg.data$group != ""),]
	group <- dat$group
	sfit <- survfit(Surv(Time, Exist) ~ group, data = dat)
	pdf(paste0(dbname,"-",gene1,"-",gene2,"-survival.pdf"),width=4,height=4.5)
	p <- ggsurvplot(sfit, pval = TRUE, risk.table = TRUE,risk.table.y.text = FALSE,risk.table.col = "strata", linetype = "solid", pval.coord = c(xsit, 1), font.x = c(12, face = "bold"),font.y = c(12, face = "bold"),  xlab = "Time (Months)", ggtheme = theme_survminer(), legend.labs=c(paste0(gene1," high + ", gene2," high"),paste0(gene1," low + ", gene2 ," low")), legend.title="", palette = c("#EE0000FF","#3B4992FF"))+guides(colour = guide_legend(nrow = 2),risk.table.fontsize=8)
	p$table <- p$table + theme(axis.line = element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
	print(p)
	dev.off()
}



#two gene survival
# CGGA IDHWT
# CGGA IDHW/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/two_genes_survival")
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/2023_revision/two_genes")
for(i in c(1:5)){
		gene1 <- genes1[i]
		gene2 <- genes2[i]
		file="/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database_CGGA/database/CGGA_mRNAseq_db_all_glioma.txt"
		mRNAseq_325_info <- read.csv("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database_CGGA/ori_data/CGGA.mRNAseq_325_clinical.20200506.xls",header=T,sep="\t")
		mRNAseq_693_info <- read.csv("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database_CGGA/ori_data/CGGA.mRNAseq_693_clinical.20200506.xls",header=T,sep="\t")
		pheno <- rbind(mRNAseq_325_info,mRNAseq_693_info)
		rownames(pheno) <- pheno$CGGA_ID
		DB <- read.csv(file,header=T,sep="\t",row.names=1)[rownames(pheno),]
		commP <- intersect(rownames(pheno),rownames(DB))
		DB <- DB[commP,]
		pheno <- pheno[commP,]
		colnames(DB) <- gsub("HLA\\.","HLA-",colnames(DB))
		pheno[which(is.na(pheno$Grade)),"Grade"] <- "no"
		pheno[which(is.na(pheno$IDH_mutation_status)),"IDH_mutation_status"] <- "no"
		Grade <- pheno[rownames(DB),"Grade"]
		IDH.status <- pheno[rownames(DB),"IDH_mutation_status"]
		DB2 <- data.frame(DB[,c(1:2)],Grade=Grade,IDH.status=IDH.status,DB[,c(3:length(colnames(DB)))])
		DB2[,"Time"] <- DB2[,"Time"]/30
		DB3 <- DB2[which(DB2$IDH.status == "Wildtype"),]
		DB <- DB3[,-c(3,4)]
		colnames(DB) <- gsub("HLA\\.","HLA-",colnames(DB))
		dbname <- "CGGA_GBM_IDHWT"
		xsit=50
		info <- DB[,c("Time","Exist",gene1,gene2)]
		deg.data <- info
		deg.data$group <- ""
		deg.data$group[which(deg.data[,gene1] >= median(info[,gene1]) & (deg.data[,gene2] >= median(info[,gene2])))] <- paste0(gene1,"-high + ",gene2,"-high")
		deg.data$group[which(deg.data[,gene1] < median(info[,gene1]) & (deg.data[,gene2] < median(info[,gene2])))] <- paste0(gene1,"-low + ",gene2,"-low")
		dat <- deg.data[which(deg.data$group != ""),]
		group <- dat$group
		sfit <- survfit(Surv(Time, Exist) ~ group, data = dat)
		pdf(paste0(dbname,"-",gene1,"-",gene2,"-survival.pdf"),width=4,height=4.5)
		p <- ggsurvplot(sfit, pval = TRUE, risk.table = TRUE,risk.table.y.text = FALSE,risk.table.col = "strata", linetype = "solid", pval.coord = c(xsit, 1), font.x = c(12, face = "bold"),font.y = c(12, face = "bold"),  xlab = "Time (Months)", ggtheme = theme_survminer(), legend.labs=c(paste0(gene1," high + ", gene2," high"),paste0(gene1," low + ", gene2 ," low")), legend.title="", palette = c("#EE0000FF","#3B4992FF"))+guides(colour = guide_legend(nrow = 2),risk.table.fontsize=8)
		p$table <- p$table + theme(axis.line = element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
		print(p)
		dev.off()
}

