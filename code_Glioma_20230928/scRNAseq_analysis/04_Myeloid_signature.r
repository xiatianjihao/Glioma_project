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
library(gridExtra)
library(ggpubr)
library(rstatix)


setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma")

Macrophage <- readRDS("RDS_file/Macrophage.rds")

levels(Macrophage) <- c("MDM-1","MDM-2","MDM-3","MG-1","MG-2","MG-3")

#signature analysis
Angiogenesis_genes <- c("CCND2","CCNE1","CD44","CXCR4","E2F3","EDN1","EZH2","FGF18","FGFR1","FYN","HEY1","ITGAV","JAG1","JAG2","MMP9","NOTCH1","PDGFA","PTK2","SPP1","STC1","TNFAIP6","TYMP","VAV2","VCAN","VEGFA")
Phagocytosis_genes <- c("SIRPB1","ABCA7","RACK1","MERTK","OLFM4","COLEC10","LMAN2","RAB31","CD300A","SPACA3","CNN2","SIRPA","CSK","CD300LF","CYBA","DNM2","DOCK2","ANO6","AHSG","F2RL1","FCER1G","FCGR2B","FCN1","FCN2","FGR","SYT11","FPR2","ALOX15","APPL1","GAS6","STAP1","GATA2","PYCARD","HCK","NCKAP1L","HMGB1","APOA1","APOA2","IFNG","IL1B","IL2RB","IL2RG","IL15","IL15RA","ITGA2","ITGAV","MYO18A","MIR17","MIR181B1","MIR183","MIR20A","MBL2","MFGE8","PLA2G5","PLCG2","PLSCR1","TREM2","APPL2","SIRPG","PIP4P2","LYAR","PRKCG","PRTN3","AZU1","CAMK1D","PTPRC","PTPRJ","PTX3","RAB27A","BCR","CCL2","NOD2","SFTPD","ATG3","CLEC7A","SLC11A1","SOD1","SYK","TGM2","TLR2","TNF","C2","C3","C4A","C4B","TUB","TULP1","COLEC11","CALR","DYSF","FCN3","SNX3","SPHK1","FER1L5","SYT7","ADIPOQ","ATG5","CD36","SCARB1","CD47")


Macrophage <- AddModuleScore(Macrophage,features=list(Angiogenesis_genes),seed=20,name="Angiogenesis_genes")
Macrophage <- AddModuleScore(Macrophage,features=list(Phagocytosis_genes),seed=20,name="Phagocytosis_genes")

df_Angiogenesis <- data.frame(group=c(rep("Angiogenesis_signature",42152)), celltype=Macrophage@meta.data[,"anno"],signature=c(Macrophage@meta.data[,"Angiogenesis_genes1"]))
compare_means(signature ~ celltype, data = df_Angiogenesis)
my_comparisons <- list(c("MDM-1","MDM-2"),c("MDM-1","MDM-3"),c("MDM-2","MDM-3"),c("MG-2","MG-1"),c("MG-1","MG-3"),c("MG-2","MG-3"))
p3 <- ggviolin(df_Angiogenesis, x = "celltype", y = "signature", color = "celltype", palette = c('#0082b8','#3db940','#a4cece',"#9f6164","#cc9f66","#fbac00"), add = c("boxplot"))+theme(title=element_text(size=12,color="black"),axis.text.x=element_text(angle=90,size=12,color="black", vjust = 0.5, hjust=1),axis.text.y=element_text(size=12,color="black"),axis.line=element_line(size=0.5,color="black"))+theme(plot.title = element_text(hjust = 0.5))+ylab("Angiogenesis Signature Score")+xlab("")+stat_compare_means(comparisons = my_comparisons,label = "p.signif")+ylim(c(-0.5,1.7))
ggsave(p3,file="Angiogenesis_signature.pdf",width=3,height=4)


df_Phagocytosis <- data.frame(group=c(rep("Phagocytosis_signature",42152)), celltype=Macrophage@meta.data[,"anno"],signature=c(Macrophage@meta.data[,"Phagocytosis_genes1"]))
compare_means(signature ~ celltype, data = df_Phagocytosis)
my_comparisons <- list(c("MDM-1","MDM-2"),c("MDM-1","MDM-3"),c("MDM-2","MDM-3"),c("MG-2","MG-1"),c("MG-1","MG-3"),c("MG-2","MG-3"))
p3 <- ggviolin(df_Phagocytosis, x = "celltype", y = "signature", color = "celltype", palette = c('#0082b8','#3db940','#a4cece',"#9f6164","#cc9f66","#fbac00"), add = c("boxplot"))+theme(title=element_text(size=12,color="black"),axis.text.x=element_text(angle=90,size=12,color="black", vjust = 0.5, hjust=1),axis.text.y=element_text(size=12,color="black"),axis.line=element_line(size=0.5,color="black"))+theme(plot.title = element_text(hjust = 0.5))+ylab("Phagocytosis Signature Score")+xlab("")+stat_compare_means(comparisons = my_comparisons,label = "p.signif")+ylim(c(-0.1,0.8))
ggsave(p3,file="Phagocytosis_signature_v2.pdf",width=3,height=4)


Macrophage <- SetIdent(Macrophage,value="anno")
levels(Macrophage) <- c("MDM-1","MDM-2","MDM-3","MG-1","MG-2","MG-3")
Coinhibitory_ICT <- unique(c("PDCD1","CD274","CTLA4","CD276","BTLA","LAG3","CEACAM1", "IDO1", "IDO2","LGALS9", "HAVCR2","TIGIT"))
breaks <- c(-1,0,1,2)
p <- DotPlot(Macrophage,features = Coinhibitory_ICT)+coord_flip()+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=12,color="black"),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 12),legend.text = element_text(size = 12))+scale_colour_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")),breaks = breaks, labels = format(breaks))
ggsave(p, file="ICT_dotplot.pdf",width=4.5,height=3.5)


##############################Difference of Grade III and Grade IV in MDM
MDM <- readRDS("RDS_file/MDM.rds")

MDM@meta.data$Grade <- "NA"
MDM@meta.data[which(MDM@meta.data$orig.ident %in% c("Glioma6","Glioma7")),]$Grade <- "GradeIII"
MDM@meta.data[which(MDM@meta.data$orig.ident %in% c("Glioma1","Glioma2","Glioma3","Glioma4","Glioma5","Glioma8","Glioma9")),]$Grade <- "GradeIV"

num_matrix=matrix(nrow=length(levels(MDM)),ncol=length(unique(MDM@meta.data$orig.ident)))
j <- 0
for(i in levels(MDM)){
    j <- j+1
    cluster=WhichCells(MDM, idents=i)
    a <- 0
    for (sample in unique(MDM@meta.data$orig.ident)){
        sample2=rownames(MDM@meta.data)[which(MDM@meta.data$orig.ident==sample)]
        num=length(intersect(cluster,sample2))
        a <- a+1
        num_matrix[j,a] <- num
    }
}
colnames(num_matrix)=unique(MDM@meta.data$orig.ident)
num_matrix2 <- rbind(num_matrix,as.vector(colSums(num_matrix)))
rownames(num_matrix2)=c(levels(MDM),"Total")
num_matrix2
num_matrix3 <- as.data.frame(num_matrix2)
num_matrix3$total <- rowSums(num_matrix2)
num_matrix3
write.table(num_matrix3, file="Number_of_eachcell_in_eachsample_split.xls",sep="\t",quote=FALSE)

percent_matrix=matrix(nrow=length(levels(MDM)),ncol=length(unique(MDM@meta.data$orig.ident)))
j <- 0
for(i in levels(MDM)){
    j <- j+1
    cluster=WhichCells(MDM, idents=i)
    a <- 0
    for (sample in unique(MDM@meta.data$orig.ident)){
        sample2=rownames(MDM@meta.data)[which(MDM@meta.data$orig.ident==sample)]
        percent=length(intersect(cluster,sample2))/length(sample2)
        a <- a+1
        percent_matrix[j,a] <- percent
    }
}
colnames(percent_matrix)=unique(MDM@meta.data$orig.ident)
rownames(percent_matrix)=c(levels(MDM))
percent_matrix
write.table(percent_matrix, file="Percent_of_eachcell_in_eachsample_split.xls",sep="\t",quote=FALSE)

file1="Percent_of_eachcell_in_eachsample_split.xls"
percent1 <- read.csv(file1,header=T,row.names=1,sep="\t")
df <- data.frame(sample = rep(colnames(percent1), each =3), celltype = rep(rownames(percent1),9), Fraction = as.numeric(as.matri
x(percent1)))
df$group <- "NA"
df$group[which(df$sample %in% c("Glioma6","Glioma7"))] <- "GradeIII"
df$group[grep("NA",df$group)] <- "GradeIV"
df$group2 <- paste0(df$celltype,"_",df$group)

stat.test <- df %>% group_by(celltype) %>% wilcox_test(Fraction ~ group) %>% adjust_pvalue(method = "bonferroni") %>% add_significance()
stat.test <- stat.test %>% add_xy_position(x = "group")

df$celltype <- factor(df$celltype,levels=c("MDM-1","MDM-2","MDM-3"))
pdf("Percent_bar_MDM_grade.pdf",width=4,height=2.5)
bxp <- ggboxplot(df, x = "group", y = "Fraction", fill = "group", palette="jco",facet.by = "celltype",nrow=1, scales = "free", add = c("jitter"), order = c("GradeIII","GradeIV"))+stat_pvalue_manual(stat.test, hide.ns = FALSE, label = "{p.adj.signif}")+scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+theme(title=element_text(size=10,color="black"),axis.text.x=element_text(angle=45,size=10,color="black",hjust=1),axis.text.y=element_text(size=10,color="black"),axis.line=element_line(size=0.3,color="black"),strip.background = element_blank(),plot.title = element_text(size = 10,vjust = 1,hjust = 0.5))+facet_wrap(~celltype, scales = "free", ncol = 7)
bxp
dev.off()


MDM <- AddModuleScore(MDM,features=list(Angiogenesis_genes),seed=20,name="Angiogenesis_genes")
MDM <- AddModuleScore(MDM,features=list(Phagocytosis_genes),seed=20,name="Phagocytosis_genes")

df_Angiogenesis <- data.frame(cells=rownames(MDM@meta.data), celltype=MDM@meta.data[,"anno"],signature=c(MDM@meta.data[,"Angiogenesis_genes1"]),Grade=MDM@meta.data[,"Grade"])

stat.test <- df_Angiogenesis %>% group_by(celltype) %>% wilcox_test(signature ~ Grade) %>% adjust_pvalue(method = "bonferroni") %>% add_significance()
stat.test <- stat.test %>% add_xy_position(x = "Grade")

df_Angiogenesis$celltype <- factor(df_Angiogenesis$celltype,levels=c("MDM-1","MDM-2","MDM-3"))
p3 <- ggviolin(df_Angiogenesis, x = "Grade", group="Grade", y = "signature", fill="Grade", palette = "jco", add = c("boxplot"),order=c("GradeIII","GradeIV"),facet.by="celltype")+stat_pvalue_manual(stat.test, hide.ns = FALSE, label = "{p.adj.signif}")+ylab("Angiogenesis Signature Score")+scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+theme(title=element_text(size=10,color="black"),axis.text.x=element_text(angle=45,size=10,color="black",hjust=1),axis.text.y=element_text(size=10,color="black"),axis.line=element_line(size=0.3,color="black"),strip.background = element_blank(),plot.title = element_text(size = 10,vjust = 1,hjust = 0.5))+facet_wrap(~celltype, scales = "free", ncol = 3)
ggsave(p3,file="MDM_Angiogenesis_signature.pdf",width=5,height=3)


df_Phagocytosis <- data.frame(cells=rownames(MDM@meta.data), celltype=MDM@meta.data[,"anno"],signature=c(MDM@meta.data[,"Phagocytosis_genes1"]),Grade=MDM@meta.data[,"Grade"])

stat.test <- df_Phagocytosis %>% group_by(celltype) %>% wilcox_test(signature ~ Grade) %>% adjust_pvalue(method = "bonferroni") %>% add_significance()
stat.test <- stat.test %>% add_xy_position(x = "Grade")

df_Phagocytosis$celltype <- factor(df_Phagocytosis$celltype,levels=c("MDM-1","MDM-2","MDM-3"))
p3 <- ggviolin(df_Phagocytosis, x = "Grade", group="Grade", y = "signature", fill="Grade", palette = "jco", add = c("boxplot"),order=c("GradeIII","GradeIV"),facet.by="celltype")+stat_pvalue_manual(stat.test, hide.ns = FALSE, label = "{p.adj.signif}")+ylab("Phagocytosis Signature Score")+scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+theme(title=element_text(size=10,color="black"),axis.text.x=element_text(angle=45,size=10,color="black",hjust=1),axis.text.y=element_text(size=10,color="black"),axis.line=element_line(size=0.3,color="black"),strip.background = element_blank(),plot.title = element_text(size = 10,vjust = 1,hjust = 0.5))+facet_wrap(~celltype, scales = "free", ncol = 3)
ggsave(p3,file="MDM_Phagocytosis_signature.pdf",width=5,height=3)




#Enrichment for macrophage subsets
library(ReactomePA)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(enrichR)
library(msigdbr)
library(tidyverse)

setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/marker_anno")

files <- list.files(path="/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/marker_anno",pattern=".xls")[2:7]

genelist <- lapply(files, function(x){
	 name <- gsub("_marker.xls","",x)
	 print(name)
	 dat <- read.csv(x,header=T,sep="\t",row.names=1)
     marker <- rownames(dat[which(dat$myAUC >0.7 & dat$avg_logFC >0),])
     change_IDs <- bitr(marker,fromType="SYMBOL",toType="ENTREZID", OrgDb="org.Hs.eg.db")
     id <- change_IDs$ENTREZID
})

genelist <- setNames(genelist, c("MDM-1","MDM-2","MDM-3","MG-1","MG-2","MG-3"))

h_gene_sets = msigdbr(species = "human", category = "H")
head(h_gene_sets)
GO_gene_sets = msigdbr(species = "human", category = "C5")
head(GO_gene_sets)
path_gene_sets = msigdbr(species = "human", category = "C2")
head(path_gene_sets)

genesets <- bind_rows(h_gene_sets,GO_gene_sets,path_gene_sets)

msigdbr_t2g = genesets %>% dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()

ck <- compareCluster(geneCluster = genelist, fun = "enricher", TERM2GENE=msigdbr_t2g)


save(ck,file="CK.rds")
data <- as.data.frame(ck)
MDMC1 <- data[which(data$Cluster == "MDM-1"),]
MDMC1_select <- MDMC1[c(1,3,5),]
MDMC2 <- data[which(data$Cluster == "MDM-2"),]
MDMC2_select <- data[c(891,893,897),]
MDMC3 <- data[which(data$Cluster == "MDM-3"),]
MDMC3_select <- data[c(1465,1466,1467),]
MGC1 <- data[which(data$Cluster == "MG-1"),]
MGC1_select <- data[c(1937,1941,1984),]
MGC2 <- data[which(data$Cluster == "MG-2"),]
MGC2_select <- data[c(2166,2172,2176),]
MGC3 <- data[which(data$Cluster == "MG-3"),]
MGC3_select <- data[c(2521,2523,2527),]


df <- rbind(MDMC1_select,MDMC2_select,MDMC3_select,MGC1_select,MGC2_select,MGC3_select)
df$GeneRatio2 <- "NA"
for(i in c(1:18)){
	df$GeneRatio2[i] <- as.numeric(strsplit(df$GeneRatio[i],"/")[[1]][1])/as.numeric(strsplit(df$GeneRatio[i],"/")[[1]][2])
}

df$Description <- factor(df$Description,levels=unique(df$Description))
df$Cluster <- factor(df$Cluster,levels=c("MDM-1","MDM-2","MDM-3","MG-1","MG-2","MG-3"))
pdf("Enrichment_result.pdf",width=11,height=5.2)
ggplot(df,aes(x = Cluster, y = Description))+geom_point(aes(size = Count, color = -log10(p.adjust)))+theme_bw(base_size = 12)+theme(text=element_text()) + scale_y_discrete(position="right")+scale_colour_gradient(low = "#FFD700", high = "red", na.value = NA)+ylab(NULL) + xlab("Cluster")+theme(axis.text = element_text(color="black"),axis.text.x = element_text(size=12,angle = 90, vjust = 0.5, hjust=1))+scale_size(range = c(2,8))
dev.off()

