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
library(ggpubr)
library(ggthemes)
library(ggrepel)

counts <- read.csv("GSE164485_gene_expression_20210104.csv",header=T,sep=",",row.names=1)
data <- CreateSeuratObject(counts = t(counts), project = "ctrl")

meta <- read.csv("GSE164485_meta_data_20210104.csv",header=T,sep=",")


data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
ribo_genes <- rownames(data@assays$RNA@data)[grep("^RP[SL]",rownames(data@assays$RNA@data))]
data[["percent.ribo"]] <- PercentageFeatureSet(data, features=ribo_genes)
hb_genes <- rownames(data@assays$RNA@data)[grep("^HB[^(P)]",rownames(data@assays$RNA@data))]
data[["percent.hb"]] <- PercentageFeatureSet(data, features=hb_genes)
HSP_genes <- rownames(data)[grep("HSP",rownames(data))]
data[["percent.hsp"]] <- PercentageFeatureSet(data, features=HSP_genes)

rownames(meta) <- meta$X
data@meta.data$subject_id2 <- meta[rownames(data@meta.data),"subject_id2"]
data@meta.data$tissue <- meta[rownames(data@meta.data),"tissue"]
data@meta.data$celltype <- meta[rownames(data@meta.data),"celltype"]

data <- SetIdent(data,value="subject_id2")

ctrl <- subset(data,idents=c("Ctrl 1","Ctrl 2","Ctrl 3","Ctrl 4"))

ctrl <- SetIdent(ctrl,value="celltype")
saveRDS(ctrl,file="ctrl.rds")

microglia <- subset(ctrl, ident="Mic")
saveRDS(microglia,file="microglia.counts.rds")
mono_maco <- subset(ctrl, ident="Mo")
saveRDS(mono_maco,file="mono_macro.counts.rds")
myeloid <- subset(ctrl,idents=c("Mic","Mo"))
saveRDS(myeloid,file="myeloid.rds")
lym <- subset(ctrl,ident="LM")
saveRDS(lym,file="Lymphocyte.counts.rds")

genes <- c("LYZ","MNDA","S100A8","CX3CR1","TMEM119","P2RY12","CD4","CD8A","CD8B")
pdf("Ctrl_vlnplot.pdf",width=12,height=6)
VlnPlot(ctrl,features=genes)
dev.off()


df <- ctrl@meta.data[,c("tissue","celltype")]

gene_expr <- t(ctrl@assays$RNA@data[genes,])

df_mat <- merge(df,gene_expr,by="row.names")

df_mat2 <- data.frame(cell=rep(df_mat[,1],9),celltype=rep(df_mat$celltype,9),Expr=c(df_mat[,"LYZ"],df_mat[,"MNDA"],df_mat[,"S100A8"],df_mat[,"CX3CR1"],df_mat[,"TMEM119"],df_mat[,"P2RY12"],df_mat[,"CD4"],df_mat[,"CD8A"],df_mat[,"CD8B"]),genes=c(rep("LYZ",26044),rep("MNDA",26044),rep("S100A8",26044),rep("CX3CR1",26044),rep("TMEM119",26044),rep("P2RY12",26044),rep("CD4",26044),rep("CD8A",26044),rep("CD8B",26044)))

pdf("Ctrl_vlnplot_ggpubr.pdf",width=12,height=8)
ggviolin(df_mat2, x="celltype", y="Expr", color="celltype", add=c("mean"),facet.by="genes")+theme(title=element_text(size=10,color="black"),axis.text.x=element_text(angle=45,size=10,color="black",hjust=1),axis.text.y=element_text(size=10,color="black"),axis.line=element_line(size=0.3,color="black"),strip.background = element_blank(),plot.title = element_text(size = 10,vjust = 1,hjust = 0.5))+facet_wrap(~genes, scales = "free", ncol = 3)
dev.off()

####################myeloid cell re-cluster#########################
data <- readRDS("myeloid.rds")
p <- VlnPlot(object = data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb","percent.hsp"), group.by="orig.ident",  ncol = 5, pt.size=0)
ggsave(p, file="Merged_sample_QC.pdf",width=40,height=5)

data<- SCTransform(data, method = "glmGamPoi", assay = "RNA",new.assay.name = "SCT", vars.to.regress = c("nCount_RNA","percent.mt"), return.only.var.genes=FALSE)
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000, assay="RNA")
saveRDS(data,file="merged_sctrans.rds")

DefaultAssay(data) <- "SCT"
length(VariableFeatures(data))

data <- RunPCA(data, features = VariableFeatures(data))
p <- PCAPlot(object = data)
ggsave(p,file="merged_PCA.pdf",width=5,height=5)

#harmony remove batch effect
p1 <- DimPlot(object = data, reduction = "pca", group.by = "subject_id2")
p2 <- VlnPlot(object = data, features = "PC_1", group.by = "subject_id2")
#plot_grid(p1,p2)
ggsave(plot_grid(p1,p2),file="All_sample_pca_before_harmony.pdf",width=14,height=7)
pdf("RunHarmony.pdf")
data <- RunHarmony(object=data, group.by.vars="subject_id2",assay.use="SCT", plot_convergence = TRUE)
dev.off()
merged_harmony_embeddings <- Embeddings(data, 'harmony')
merged_harmony_embeddings[1:5, 1:5]
p1 <- DimPlot(object = data, reduction = "harmony", group.by = "subject_id2")
p2 <- VlnPlot(object = data, features = "harmony_1", group.by = "subject_id2")
#plot_grid(p1,p2)
ggsave(plot_grid(p1,p2),file="All_sample_pca_after_harmony_group_by_sample.pdf",width=16,height=8)
p <- ElbowPlot(data,reduction="harmony",ndims=50)
ggsave(p, file="ElbowPlot.pdf",width=5,height=5)

data <- FindNeighbors(data, dims = 1:30, reduction = "harmony", verbose = FALSE)
data <- FindClusters(data, reduction.type = "harmony",resolution = seq(0.1, 1.0, 0.1),dims.use = 1:30, verbose = FALSE)
data <- RunUMAP(data, dims = 1:30,reduction = "harmony",verbose = FALSE)
data <- RunTSNE(data, dims = 1:30,reduction = "harmony",verbose = FALSE)
saveRDS(data,file="ctrl_merged.rds")


head(data@meta.data)
sapply(grep("^SCT_snn_res",colnames(data@meta.data),value = TRUE),
       function(x) length(unique(data@meta.data[,x])))

data <- SetIdent(data, value="SCT_snn_res.1")
p1 <- DimPlot(data, reduction = "umap", group.by = "subject_id2")
p2 <- DimPlot(object = data, reduction = 'umap', label=TRUE, label.size=5)
ggsave(plot_grid(p1,p2),file="merged_umap.pdf",width=14,height=6)
p <- DimPlot(object = data, reduction = 'umap', label=TRUE, label.size=5,split.by="subject_id2",ncol=10)
ggsave(p, file="merged_umap_split.pdf",width=40,height=5)
p1 <- DimPlot(data, reduction = "tsne", group.by = "subject_id2")
p2 <- DimPlot(object = data, reduction = 'tsne', label=TRUE, label.size=5)
ggsave(plot_grid(p1,p2),file="merged_tsne.pdf",width=14,height=6)
p <- DimPlot(object = data, reduction = 'tsne', label=TRUE, label.size=5,split.by="subject_id2",ncol=10)
ggsave(p, file="merged_tsne_split.pdf",width=40,height=5)

genes <- c("LYZ","MNDA","S100A9","S100A8","S100A6","FCN1","VCAN","ITGA4","MRC1","TGFBI","THBD","CX3CR1","TMEM119","ITGAX","P2RY12","BIN1","NAV3","IFITM3","ISG15","TNF","LY6E","HLA-DMA","HLA-DMB","HLA-DQA2","HFE","MIF","ENO1","ALDOA","LDHA","EPAS1","CCL4","CCL4L2","CCL3","CCL3L1","SPP1","HSPA6","HSPA1A","HSPA1B","DNAJB1","IFI27","IFI6","MX1","IFI44L")
breaks <- c(-1.5,0,1.5)
p <- DotPlot(data,features = genes)+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=12,color="black"),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 12),legend.text = element_text(size = 12))+scale_colour_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")),breaks = breaks, labels = format(breaks))+scale_size(range = c(4,8))
ggsave(p, file="marker_dotplot_macro.pdf",width=12,height=3)

for (i in levels(data)){
  marker <- FindMarkers(data,ident.1=i, ident.2=NULL, test.use="roc",only.pos=TRUE)
  write.table(marker,file=paste0(i,"_marker.xls"),sep="\t",quote=FALSE)
}


data <- SetIdent(data, value="SCT_snn_res.1")
data@meta.data$anno <- "NA"
data@meta.data[which(data@meta.data$SCT_snn_res.1 %in% c(0,7,9)),]$anno <- "Monocyte"
data@meta.data[which(data@meta.data$SCT_snn_res.1 %in% c(1,4)),]$anno <- "MG-1"
data@meta.data[which(data@meta.data$SCT_snn_res.1 %in% c(5,6)),]$anno <- "MG-2"
data@meta.data[which(data@meta.data$SCT_snn_res.1 %in% c(2,3,8,10,11)),]$anno <- "MG-3"

data@meta.data$celltype2 <- "NA"
data@meta.data[which(data@meta.data$SCT_snn_res.1 %in% c(0,7,9)),]$celltype2 <- "Monocyte"
data@meta.data[which(data@meta.data$SCT_snn_res.1 %in% c(1,2,3,4,5,6,8,10,11)),]$celltype2 <- "MG"
saveRDS(data,file="ctrl.anno.rds")

data <- SetIdent(data,value="celltype2")
genes <- c("LYZ","FCN1","S100A8","S100A9","CX3CR1","TMEM119","P2RY12","BIN1","NAV3")
breaks <- c(-1,-0.5,0,0.5,1)
p <- DotPlot(data,features = genes)+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=12,color="black"),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 12),legend.text = element_text(size = 12))+scale_colour_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")),breaks = breaks, labels = format(breaks))+scale_size(range = c(4,8))
ggsave(p, file="marker_dotplot_celltype.pdf",width=5,height=2)


#####################Differentially expressed genes in myeloid cells between control and tumor
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/myeloid_gbm_vs_ctrl_2023")

Glioma <- readRDS("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/Myeloid_anno.rds")

Ctrl <- readRDS("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/myeloid_gbm_vs_ctrl_2023/ctrl.anno.rds")
rownames(Ctrl@assays$RNA@counts) <- gsub("\\.","-",rownames(Ctrl@assays$RNA@counts))
rownames(Ctrl@assays$RNA@data) <- gsub("\\.","-",rownames(Ctrl@assays$RNA@data))
rownames(Ctrl@assays$RNA@scale.data) <- gsub("\\.","-",rownames(Ctrl@assays$RNA@scale.data))
Glioma@meta.data$group <- "Glioma"
Ctrl@meta.data$group <- "Ctrl"

data <- merge(Glioma, y = Ctrl, add.cell.ids = c("Glioma","Ctrl"), project = "merged")
data <- SetIdent(data,value="celltype2")
saveRDS(data,file="Glioma_vs_Ctrl.rds")

options(future.globals.maxSize = 8000 * 1024^2)
for(i in c("Monocyte","MG")){
    cluster <- subset(data,idents=i)
    Idents(cluster) <- "group"
    diff_expr <- FindMarkers(cluster,ident.1="Ctrl",ident.2="Glioma",assay="RNA",test.use = "wilcox",logfc.threshold=0,verbose=FALSE)
    write.table(diff_expr,file=paste0(i,"_diff_expr.xls"),sep="\t",quote=FALSE)
}

#Enrichment for macrophage subsets
library(ReactomePA)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(enrichR)
library(msigdbr)
library(tidyverse)

setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/control_brain_2023/GSE164485/control_brain_analysis/myeloid")

files <- list.files(path="/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/control_brain_2023/GSE164485/control_brain_analysis/myeloid",pattern=".xls")

genelist <- lapply(files, function(x){
    name <- gsub("_diff_expr.xls","",x)
    print(name)
    dat <- read.csv(x,header=T,sep="\t",row.names=1)
     marker <- rownames(dat[which(dat$avg_log2FC >1),])
     change_IDs <- bitr(marker,fromType="SYMBOL",toType="ENTREZID", OrgDb="org.Hs.eg.db")
     id <- change_IDs$ENTREZID
})

genelist <- setNames(genelist, c("MG","Monocyte"))

h_gene_sets = msigdbr(species = "human", category = "H")
head(h_gene_sets)
GO_gene_sets = msigdbr(species = "human", category = "C5")
head(GO_gene_sets)

genesets <- bind_rows(h_gene_sets,GO_gene_sets)

msigdbr_t2g = genesets %>% dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()

ck <- compareCluster(geneCluster = genelist, fun = "enricher", TERM2GENE=msigdbr_t2g)


save(ck,file="CK.Glioma.up.rds")
data <- as.data.frame(ck)
Monocyte <- data[which(data$Cluster == "Monocyte"),]
Monocyte_select <- Monocyte[c(1:6,9),]
MGC <- data[which(data$Cluster == "MG"),]
MGC_select <- MGC[c(1:5,8:9),]

Monocyte_select$GeneRatio2 <- "NA"
for(i in c(1:7)){
        Monocyte_select$GeneRatio2[i] <- as.numeric(strsplit(Monocyte_select$GeneRatio[i],"/")[[1]][1])/as.numeric(strsplit(Monocyte_select$GeneRatio[i],"/")[[1]][2])
}

Monocyte_select <- Monocyte_select[rev(order(Monocyte_select$p.adjust)),]
Monocyte_select$Description <- factor(Monocyte_select$Description,levels=unique(Monocyte_select$Description))
pdf("Monocyte_Glioma_up_Enrichment_barplot.pdf",width=7,height=2.5)
ggplot(Monocyte_select,aes(x = -log10(p.adjust), y = Description))+geom_bar(stat="identity",color="orange",fill="orange")+theme_classic(base_size = 12)+theme(text=element_text())+ylab(NULL) + xlab("-log10(p.adjust)")+theme(axis.text = element_text(color="black"),axis.text.x = element_text(size=12,angle = 90, vjust = 0.5, hjust=1))+scale_size(range = c(2,8))
dev.off()

MGC_select$GeneRatio2 <- "NA"
for(i in c(1:7)){
        MGC_select$GeneRatio2[i] <- as.numeric(strsplit(MGC_select$GeneRatio[i],"/")[[1]][1])/as.numeric(strsplit(MGC_select$GeneRatio[i],"/")[[1]][2])
}

MGC_select <- MGC_select[rev(order(MGC_select$p.adjust)),]
MGC_select$Description <- factor(MGC_select$Description,levels=unique(MGC_select$Description))
pdf("MG_Glioma_up_Enrichment_barplot.pdf",width=7,height=2.5)
ggplot(MGC_select,aes(x = -log10(p.adjust), y = Description))+geom_bar(stat="identity",color="orange",fill="orange")+theme_classic(base_size = 12)+theme(text=element_text())+ylab(NULL) + xlab("-log10(p.adjust)")+theme(axis.text = element_text(color="black"),axis.text.x = element_text(size=12,angle = 90, vjust = 0.5, hjust=1))+scale_size(range = c(2,8))
dev.off()

