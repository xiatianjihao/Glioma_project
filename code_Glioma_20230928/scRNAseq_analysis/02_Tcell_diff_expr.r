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

setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Tcell/Glioma_vs_PBMC_diff")

Glioma <- readRDS("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Tcell/Glioma/GBM_Tcell_anno.rds")

PBMC <- readRDS("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Tcell/PBMC/PBMC_Tcell_anno.rds")
Glioma@meta.data$group <- "Glioma"
PBMC@meta.data$group <- "PBMC"


data <- merge(Glioma, y = PBMC, add.cell.ids = c("Glioma","PBMC"), project = "merged")
data <- SetIdent(data,value="celltype")
saveRDS(data,file="Glioma_vs_PBMC_Tcell.rds")

options(future.globals.maxSize = 8000 * 1024^2)
for(i in c("CD4_Naive","CD4_Stress","Treg","CD8_Cytotoxicity")){
    cluster <- subset(data,idents=i)
    Idents(cluster) <- "group"
    diff_expr <- FindMarkers(cluster,ident.1="Glioma",ident.2="PBMC",assay="RNA",test.use = "wilcox",logfc.threshold=0,verbose=FALSE)
    write.table(diff_expr,file=paste0(i,"_diff_expr.xls"),sep="\t",quote=FALSE)
}

#Volcano plot
files <- list.files(path="/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Tcell/Glioma_vs_PBMC_diff",pattern="*.xls")

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
	if(name == "CD4_Naive"){
		up.genes <- c("CREM","FOS","JUN","GZMK","LDHA","HSPA1A")
		down.genes <- c("TCF7","TXNIP")
	}else if(name == "CD4_Stress"){
		up.genes <- c("HSPA1A","HSPA1B","CAPG","IFNG","CXCR4","CTSL","COTL1")
		down.genes <- c("KLF2","AES","IER2","TXNIP")
	}else if(name == "Treg"){
		up.genes <- c("TNFRSF9","DUSP2","TNFRSF18","IL2RA", "IL2RB","IL2RG","ICOS","DUSP4","CDKN1A","CREM","CTLA4")
		down.genes <- c("MAL","KLF2","RPL34","KLF3","RPS27","AES","GIMAP7")
	}else{
		up.genes <- c("NR4A1","NR4A2","DUSP2", "PDE4B", "SAT1", "CXCR4", "HSPA1A", "DUSP1","TSC22D3")
		down.genes <- c("AES","TXNIP","GIMAP7","SEPT7","PLEK","PRF1","GZMH")
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

#signature analysis
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Tcell/Glioma_vs_PBMC_diff")

data <- readRDS("Glioma_vs_PBMC_Tcell.rds")
Cytotoxicity_markers <- c("NKG7","GNLY","PRF1","GZMA","GZMB","GZMH")
Exhaustion_markers <- c("PDCD1","TIGIT","LAG3","HAVCR2","CTLA4","LAYN","TOX","TOX2","ENTPD1","ITGAE","IFNG","CXCL13")

data <- AddModuleScore(data, list(Cytotoxicity_markers), name = "Cytotoxicity", seed = 2)
data <- AddModuleScore(data, list(Exhaustion_markers), name = "Exhaustion", seed = 2)
Cytotoxicity_df <- data@meta.data[,c("anno","group","Cytotoxicity1")]
colnames(Cytotoxicity_df) <- c("CellType","Group","Signature score")
Cytotoxicity_df$CellType <- factor(Cytotoxicity_df$CellType,levels=c("CD4_Naive","CD4_Memory","CD4_Effector_memory","CD4_Stress","Treg","CD8_Naive","CD8_Cytotoxicity","CD8_Effector_memory","CD8_Exhaustion","CD8_Stress"))

pdf("Cytotoxicity_signature_score.pdf",width=5,height=3.5)
ggviolin(Cytotoxicity_df, x = "CellType", y = "Signature score", color = "Group", palette=c("#D73027","#4575B4"),
         add = "boxplot", add.params = list(fill = "white"))+ylab("Cytotoxicity Signature Score")+xlab("")+theme(title=element_text(size=10,color="black"),axis.text.x=element_text(angle=30,size=10,color="black", vjust = 1, hjust=1),axis.text.y=element_text(size=10,color="black"),axis.line=element_line(size=0.5,color="black"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5),plot.title = element_text(hjust = 0.5))+stat_compare_means(aes(group=Group), label = "p.signif",label.y=2.5)+ylim(c(-1.2,2.8))
dev.off()


Exhaustion_df <- data@meta.data[,c("anno","group","Exhaustion1")]
colnames(Exhaustion_df) <- c("CellType","Group","Signature score")
Exhaustion_df$CellType <- factor(Exhaustion_df$CellType,levels=c("CD4_Naive","CD4_Memory","CD4_Effector_memory","CD4_Stress","Treg","CD8_Naive","CD8_Cytotoxicity","CD8_Effector_memory","CD8_Exhaustion","CD8_Stress"))

pdf("Exhaustion_signature_score.pdf",width=5,height=3.5)
ggviolin(Exhaustion_df, x = "CellType", y = "Signature score", color = "Group", palette=c("#D73027","#4575B4"),
         add = "boxplot", add.params = list(fill = "white"))+ylab("Exhaustion Signature Score")+xlab("")+theme(title=element_text(size=10,color="black"),axis.text.x=element_text(angle=30,size=10,color="black", vjust = 1, hjust=1),axis.text.y=element_text(size=10,color="black"),axis.line=element_line(size=0.5,color="black"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5),plot.title = element_text(hjust = 0.5))+stat_compare_means(aes(group=Group), label = "p.signif",label.y=0.8)+ylim(c(-0.2,0.9))
dev.off()

#AHR differential expression
data <- SetIdent(data,value="anno")
data2 <- subset(data,ident="CD8_Cytotoxicity")
genes <- c("AHR")

df <- data2@meta.data[,c("anno","group")]
df$Expression <- data2@assays$RNA@data["AHR",]
df$Gene <- "AHR"
colnames(df) <- c("CellType","Group","Expression","Gene")

df$Group <- gsub("GBM","Glioma",df$Group)

pdf("AHR_Expression.pdf",width=2.8,height=2.8)
ggviolin(df, x = "Group", y = "Expression", color = "Group", palette=c("#D73027","#4575B4"),
         add = "mean_se")+ylab("AHR Expression")+xlab("")+theme(title=element_text(size=12,color="black"),axis.text.x=element_text(size=12,color="black", vjust = 1, hjust=1),axis.text.y=element_text(size=12,color="black"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1),plot.title = element_text(hjust = 0.5))+stat_compare_means(aes(group=Group), label = "p.signif")
dev.off()
