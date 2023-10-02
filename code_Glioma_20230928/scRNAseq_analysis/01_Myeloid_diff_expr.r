library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(Matrix)
library(cowplot)
library(dplyr)
library(harmony)
library(reshape2)
library(ggsci)
library(pheatmap)
library(ggpubr)
library(ggthemes)
library(ggrepel)



setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma_vs_PBMC_diff")

Glioma <- readRDS("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/Myeloid_anno.rds")

PBMC <- readRDS("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/PBMC/PBMC_myeloid_anno.rds")
Glioma@meta.data$group <- "Glioma"
PBMC@meta.data$group <- "PBMC"


data <- merge(Glioma, y = PBMC, add.cell.ids = c("Glioma","PBMC"), project = "merged")
data <- SetIdent(data,value="celltype")
saveRDS(data,file="Glioma_vs_PBMC.rds")

options(future.globals.maxSize = 8000 * 1024^2)
for(i in c("Monocyte","DC","MDM")){
    cluster <- subset(data,idents=i)
    Idents(cluster) <- "group"
    diff_expr <- FindMarkers(cluster,ident.1="Glioma",ident.2="PBMC",assay="RNA",test.use = "wilcox",logfc.threshold=0,verbose=FALSE)
    write.table(diff_expr,file=paste0(i,"_diff_expr.xls"),sep="\t",quote=FALSE)
}

#Volcano plot
files <- list.files(path="/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma_vs_PBMC_diff",pattern="*.xls")

for(file in files){
	name <- gsub("_diff_expr.xls","",file)
	file2 <- paste0(name,"_diff.pdf")
	res <- read.csv(file,header=T,sep="\t",row.names=1)
	res$log2fc <- log2(exp(res$avg_logFC))
	write.table(res,file=paste0(name,"_diff_expr2.xls"),sep="\t",quote=FALSE)
	data <- res
	data$logP <- -log10(data$p_val_adj)
	data$group <- "Not-significant"
	data$group[which(data$p_val_adj < 0.05 & data$log2fc >1)] <- "Up-regulated"
	data$group[which(data$p_val_adj < 0.05 & data$log2fc < -1)] <- "Down-regulated"
	data$color <- "#BBBBBB"
	data$color[which(data$p_val_adj < 0.05 & data$log2fc >1)] <- "#CC0000"
	data$color[which(data$p_val_adj < 0.05 & data$log2fc < -1)] <- "#2f5688"
	data$Label <- ""
	data <- data[order(data$log2fc, decreasing = TRUE),]
	if(name == "Monocyte"){
		up.genes <- c("HIF1A","LDHA","AREG","EREG","HBEGF","CDKN1A","IL1B","CCL4")
		down.genes <- c("CYP1B1","IGKC","LRRK2","MARCH1","POU2F2","RIPOR2")
	}else if(name == "DC"){
		up.genes <- c("CCL4","CCL3","HIF1A","HSP90AA1","YWHAH","CTSD","YWHAH","HSPB1")
		down.genes <- c("RPL34","RPL39","LGALS2","POU2F2")
	}else{
		up.genes <- c("CCL2", "CXCL3","IL1B","NUPR1","MIF","HIF1A","HLA-DPB1")
		down.genes <- c("CD52","RPL34","CD48","CLEC12A","LYZ","CSTA","FCGR3A")
	}

	top10_genes <- c(as.character(up.genes),as.character(down.genes))
	data$Label[match(top10_genes,rownames(data))] <- top10_genes
	p <- ggscatter(data, x="log2fc", y="logP", color="group", palette=c("#2f5688","#BBBBBB","#CC0000"),size=0.5,font.label=8,repel=T, xlab="log2FoldChange",ylab="-log10(Adjust P-value)")+theme_base()+geom_hline(yintercept=1.30, linetype="dashed")+geom_vline(xintercept=c(-1,1),linetype="dashed")+geom_text_repel(label=data$Label, size=3, colour = factor(data$color), box.padding = 0, max.overlaps = Inf,segment.size  = 0.1,segment.color = "black")
	ggsave(p,file=paste0(name,"_diff_expr.pdf"),width=6.5,height=4.5)
}


#generate gesa rnk file
for (file in files){
	data <- read.csv(file,header=T,sep="\t", row.names=1)
	name <- gsub("_diff_expr.xls","",file)
	data <- data[order(data$avg_logFC,decreasing = TRUE),]
	lfc <- as.data.frame(data$avg_logFC)
	rownames(lfc) <- rownames(data)
	colnames(lfc) <- "#ignore"
	write.table(lfc,file=paste0(name,".rnk"),sep="\t",quote=FALSE)
}

