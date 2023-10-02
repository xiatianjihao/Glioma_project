#First calculate signature score using ssGSVA function in GSVA package, then conduct survival analysis using signature score

setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/marker_anno")

library(GSVA)
library(GSEABase)
library(RColorBrewer)
library(pheatmap)
library(survival)
library(survminer)
library(grid)
library(gridExtra)

#TCGA GBMLGG
file <- "/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database/download_GlioVis/2023-07-12_TCGA_GBMLGG_expression.txt"
pheno_file <- "/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database/download_GlioVis/2023-07-12_TCGA_GBMLGG_pheno.txt"
pheno <- read.csv(pheno_file,header=T,sep="\t",row.names=1)
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
DB_mat <- t(DB3[,-c(1,2,3,4)])
gene_name <- rownames(DB_mat)
DB_mat2 <- data.frame(gene_name,DB_mat)
        #geneset file
geneSets <- getGmt("MDM_vs_MG.gmt")

name=DB_mat2[,1]
index <- duplicated(DB_mat2[,1])
fildup=DB_mat2[!index,]
exp=fildup[,-1]
row.names(exp)=name
        #将数据框转换成矩阵
mydata= as.matrix(exp)
res_es <- gsva(mydata, geneSets, min.sz=1, max.sz=1000, verbose=FALSE, parallel.sz=1)

colnames(res_es) <- gsub("\\.","-",colnames(res_es))

res <- t(res_es)

sur_data <- DB3[,c(1,2)]
rownames(sur_data) <- gsub("\\.","-",rownames(sur_data))
sur_data2 <- merge(sur_data, res, by="row.names")

genelist <- rownames(res_es)

dbname <- "TCGA_GBM_LGG_IDHWT"
xsit=75
for(gene in genelist){
    if(gene %in% colnames(sur_data2)){
        site=which(colnames(sur_data2) == gene)
        deg.data <- sur_data2
        deg.data$group <- ""

        deg.data$group[which(deg.data[,gene] >= median(deg.data[,gene]))] <- paste0(gene,"signature high")
        deg.data$group[which(deg.data[,gene] < median(deg.data[,gene]))] <- paste0(gene,"signature low")

        dat <- deg.data[which(deg.data$group != ""),]
        group <- dat$group
        sfit <- survfit(Surv(Time, Exist) ~ group, data = dat)
        pdf(paste0(dbname,"-",gene,"-survival.pdf"),width=4,height=4.5)
        p <- ggsurvplot(sfit, pval = TRUE, risk.table = TRUE,risk.table.y.text = FALSE,risk.table.col = "strata", linetype = "solid", pval.coord = c(xsit, 1), font.x = c(12, face = "bold"),font.y = c(12, face = "bold"),  xlab = "Time (Months)", ggtheme = theme_survminer(),legend.labs=c(paste0(gene," signature high"),paste0(gene," signature low")), legend.title="", palette = c("#EE0000FF","#3B4992FF"))+guides(colour = guide_legend(nrow = 2))
        p$table <- p$table + theme(axis.line = element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
        print(p)
        dev.off()
    }
}


##### CGGA

#IDH WT
file="/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database_CGGA/database/CGGA_mRNAseq_db_all_glioma.txt"

mRNAseq_325_info <- read.csv("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database_CGGA/ori_data/CGGA.mRNAseq_325_clinical.20200506.xls",header=T,sep="\t")
mRNAseq_693_info <- read.csv("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/GBM_survival/database_CGGA/ori_data/CGGA.mRNAseq_693_clinical.20200506.xls",header=T,sep="\t")
pheno <- rbind(mRNAseq_325_info,mRNAseq_693_info)
rownames(pheno) <- pheno$CGGA_ID

DB <- read.csv(file,header=T,sep="\t",row.names=1)
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
DB_mat <- t(DB3[,-c(1,2,3,4)])
gene_name <- rownames(DB_mat)
DB_mat2 <- data.frame(gene_name,DB_mat)
        #geneset file
geneSets <- getGmt("MDM_vs_MG.gmt")

name=DB_mat2[,1]
index <- duplicated(DB_mat2[,1])
fildup=DB_mat2[!index,]
exp=fildup[,-1]
row.names(exp)=name
        #将数据框转换成矩阵
mydata= as.matrix(exp)
res_es <- gsva(mydata, geneSets, min.sz=1, max.sz=1000, verbose=FALSE, parallel.sz=1)

colnames(res_es) <- gsub("\\.","-",colnames(res_es))

res <- t(res_es)

sur_data <- DB3[,c(1,2,3)]
rownames(sur_data) <- gsub("\\.","-",rownames(sur_data))
sur_data2 <- merge(sur_data, res, by="row.names")

genelist <- rownames(res_es)

dbname <- "CGGA_GBM_IDHWT"
xsit=50
for(gene in genelist){
    if(gene %in% colnames(sur_data2)){
        site=which(colnames(sur_data2) == gene)
        deg.data <- sur_data2
        deg.data$group <- ""
        deg.data$group[which(deg.data[,gene] >= median(deg.data[,gene]))] <- paste0(gene,"signature high")
         deg.data$group[which(deg.data[,gene] < median(deg.data[,gene]))] <- paste0(gene,"signature low")

        dat <- deg.data[which(deg.data$group != ""),]
        group <- dat$group
        sfit <- survfit(Surv(Time, Exist) ~ group, data = dat)
        pdf(paste0(dbname,"-",gene,"-survival.pdf"),width=4,height=4.5)
        p <- ggsurvplot(sfit, pval = TRUE, risk.table = TRUE,risk.table.y.text = FALSE,risk.table.col = "strata", linetype = "solid", pval.coord = c(xsit, 1), font.x = c(12, face = "bold"),font.y = c(12, face = "bold"),  xlab = "Time (Months)", ggtheme = theme_survminer(),legend.labs=c(paste0(gene," signature high"),paste0(gene," signature low")), legend.title="", palette = c("#EE0000FF","#3B4992FF"))+guides(colour = guide_legend(nrow = 2))
        p$table <- p$table + theme(axis.line = element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
        print(p)
        dev.off()
    }
}




