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

setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Tcell/PBMC")

Tcell <- readRDS("../../RDS_file/T.rds")
Tcell <- SetIdent(Tcell,value="group")
PBMC <- subset(Tcell,ident="Glioma")

counts <- PBMC@assays$RNA@counts
data <- CreateSeuratObject(counts = counts, project = "GBM", min.cells = 3, min.features = 200)

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
ribo_genes <- rownames(data@assays$RNA@data)[grep("^RP[SL]",rownames(data@assays$RNA@data))]
data[["percent.ribo"]] <- PercentageFeatureSet(data, features=ribo_genes)
hb_genes <- rownames(data@assays$RNA@data)[grep("^HB[^(P)]",rownames(data@assays$RNA@data))]
data[["percent.hb"]] <- PercentageFeatureSet(data, features=hb_genes)

p1 <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo","percent.hb"), ncol = 5)
ggsave(p1,file="Tcell_QC.pdf",width=25,height=8)


data <- SCTransform(data, method = "glmGamPoi", assay = "RNA",new.assay.name = "SCT", vars.to.regress = c("nCount_RNA","percent.mt"), return.only.var.genes=FALSE)
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000, assay="RNA")
saveRDS(data,file="Tcell_sctrans.rds")

data <- RunPCA(data, features = VariableFeatures(object = data))
#harmony remove batch effect
p1 <- DimPlot(object = data, reduction = "pca", group.by = "orig.ident")
p2 <- VlnPlot(object = data, features = "PC_1", group.by = "orig.ident")
ggsave(plot_grid(p1,p2),file="All_sample_pca_before_harmony.pdf",width=14,height=7)
pdf("harmony.pdf")
data <- RunHarmony(object=data, group.by.vars="orig.ident",assay.use="SCT", plot_convergence = TRUE)
dev.off()
data_harmony_embeddings <- Embeddings(data, 'harmony')
data_harmony_embeddings[1:5, 1:5]
p1 <- DimPlot(object = data, reduction = "harmony", group.by = "orig.ident")
p2 <- VlnPlot(object = data, features = "harmony_1", group.by = "orig.ident")
ggsave(plot_grid(p1,p2),file="All_sample_pca_after_harmony_group_by_sample.pdf",width=16,height=8)
p <- ElbowPlot(data,reduction="harmony",ndims=50)
ggsave(p, file="ElbowPlot.pdf",width=5,height=5)


#20 ok
data <- FindNeighbors(data, dims = 1:30, reduction = "harmony", verbose = FALSE)
data <- FindClusters(data, reduction.type = "harmony",resolution = seq(0.1, 1, 0.1),dims.use = 1:30, verbose = FALSE)
data <- RunUMAP(data, dims = 1:30,reduction = "harmony",  verbose = FALSE)
data <- RunTSNE(data, dims = 1:30,reduction = "harmony")
saveRDS(data,file="Tcell_with_coordinate_PC30.rds")

head(data@meta.data)
sapply(grep("^SCT_snn_res",colnames(data@meta.data),value = TRUE),
       function(x) length(unique(data@meta.data[,x])))

data <- SetIdent(data, value="SCT_snn_res.1")
p1 <- DimPlot(data, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object = data, reduction = 'umap', label=TRUE, label.size=5)
ggsave(plot_grid(p1,p2),file="data_umap.pdf",width=14,height=6)
p <- DimPlot(object = data, reduction = 'umap', label=TRUE, label.size=5,split.by="orig.ident",ncol=10)
ggsave(p, file="data_umap_split.pdf",width=40,height=5)
p1 <- DimPlot(data, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(object = data, reduction = 'tsne', label=TRUE, label.size=5)
ggsave(plot_grid(p1,p2),file="data_tsne.pdf",width=14,height=6)
p <- DimPlot(object = data, reduction = 'tsne', label=TRUE, label.size=5,split.by="orig.ident",ncol=10)
ggsave(p, file="data_tsne_split.pdf",width=40,height=5)



genes <- c("CD3D","CD3E","CD3G","CD8A","CD8B","CD4","CD40LG","CST3","CD68","AIF1","GZMA","KLRD1","NKG7","KLRB1","KLRC1")
p <- DoHeatmap(data,features = genes,size=5,draw.lines=TRUE,combin=TRUE)+scale_fill_gradientn(colors = c("blue", "white", "red"),na.value="white")+theme(text = element_text(size=11),axis.text = element_text(color="black"))
ggsave(p, file="marker_heatmap.pdf",width=20,height=5)
genes <- c("CD3D","CD3E","CD3G","CD8A","CD8B","CD4","CD40LG","CST3","CD68","AIF1","GZMA","KLRD1","NKG7","KLRB1","KLRC1")
p <- DotPlot(data,features = genes)+coord_flip()+theme(text = element_text(size=11),axis.text = element_text(color="black"))
ggsave(p, file="marker_dotplot.pdf",width=20,height=5)


#remove 1,8,9,18,21,22
data1 <- subset(data,idents=c(0,2:7,10:17,19:20))
counts <- data1@assays$RNA@counts

data <- CreateSeuratObject(counts = counts, project = "GBM", min.cells = 3, min.features = 200)

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
ribo_genes <- rownames(data@assays$RNA@data)[grep("^RP[SL]",rownames(data@assays$RNA@data))]
data[["percent.ribo"]] <- PercentageFeatureSet(data, features=ribo_genes)
hb_genes <- rownames(data@assays$RNA@data)[grep("^HB[^(P)]",rownames(data@assays$RNA@data))]
data[["percent.hb"]] <- PercentageFeatureSet(data, features=hb_genes)

p1 <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo","percent.hb"), ncol = 5)
ggsave(p1,file="Tcell_QC.pdf",width=25,height=8)

data <- subset(data, subset = percent.mt < 12)

data <- SCTransform(data, method = "glmGamPoi", assay = "RNA",new.assay.name = "SCT", vars.to.regress = c("nCount_RNA","percent.mt"), return.only.var.genes=FALSE)
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000, assay="RNA")
saveRDS(data,file="Tcell_sctrans_v2.rds")


data <- RunPCA(data, features = VariableFeatures(object = data))
#harmony remove batch effect
p1 <- DimPlot(object = data, reduction = "pca", group.by = "orig.ident")
p2 <- VlnPlot(object = data, features = "PC_1", group.by = "orig.ident")
ggsave(plot_grid(p1,p2),file="All_sample_pca_before_harmony.pdf",width=14,height=7)
pdf("harmony.pdf")
data <- RunHarmony(object=data, group.by.vars="orig.ident",assay.use="SCT", plot_convergence = TRUE)
dev.off()
data_harmony_embeddings <- Embeddings(data, 'harmony')
data_harmony_embeddings[1:5, 1:5]
p1 <- DimPlot(object = data, reduction = "harmony", group.by = "orig.ident")
p2 <- VlnPlot(object = data, features = "harmony_1", group.by = "orig.ident")
ggsave(plot_grid(p1,p2),file="All_sample_pca_after_harmony_group_by_sample.pdf",width=16,height=8)
p <- ElbowPlot(data,reduction="harmony",ndims=50)
ggsave(p, file="ElbowPlot.pdf",width=5,height=5)

data <- FindNeighbors(data, dims = 1:30, reduction = "harmony", verbose = FALSE)
data <- FindClusters(data, reduction.type = "harmony",resolution = seq(0.1, 1, 0.1),dims.use = 1:30, verbose = FALSE)
data <- RunUMAP(data, dims = 1:30,reduction = "harmony",  verbose = FALSE)
data <- RunTSNE(data, dims = 1:30,reduction = "harmony")
saveRDS(data,file="Tcell_with_coordinate_PC30_v2.rds")

head(data@meta.data)
sapply(grep("^SCT_snn_res",colnames(data@meta.data),value = TRUE),
       function(x) length(unique(data@meta.data[,x])))

data <- SetIdent(data, value="SCT_snn_res.1")
p1 <- DimPlot(data, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object = data, reduction = 'umap', label=TRUE, label.size=5)
ggsave(plot_grid(p1,p2),file="data_umap.pdf",width=14,height=6)
p <- DimPlot(object = data, reduction = 'umap', label=TRUE, label.size=5,split.by="orig.ident",ncol=10)
ggsave(p, file="data_umap_split.pdf",width=40,height=5)
p1 <- DimPlot(data, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(object = data, reduction = 'tsne', label=TRUE, label.size=5)
ggsave(plot_grid(p1,p2),file="data_tsne.pdf",width=14,height=6)
p <- DimPlot(object = data, reduction = 'tsne', label=TRUE, label.size=5,split.by="orig.ident",ncol=10)
ggsave(p, file="data_tsne_split.pdf",width=40,height=5)

for (i in levels(data)){
  marker <- FindMarkers(data,ident.1=i, ident.2=NULL, test.use="roc",only.pos=TRUE)
  write.table(marker,file=paste0(i,"_marker.xls"),sep="\t",quote=FALSE)
}

genes <- c("CD3D","CD3E","CD3G","CD8A","CD8B","CD4","CD40LG","CST3","CD68","AIF1","GZMA","KLRD1","NKG7","KLRB1","KLRC1")
breaks <- c(-1,0,1,2)
p <- DotPlot(data,features = genes)+coord_flip()+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=12,color="black"),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 12),legend.text = element_text(size = 12))+scale_colour_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")),breaks = breaks, labels = format(breaks))+ scale_size(range = c(4, 8))
ggsave(p, file="marker_dotplot.pdf",width=13,height=5)


heat_genes <- c("CD4","TCF7","CCR7","SELL","LEF1","CD28","ITGB1","S100A4","LGALS1","DNAJA1","DNAJB1","HSP90AB1","HSP90AA1","FOXP3","TIGIT","CTLA4","CD8A","CD8B","PRF1","GZMB","GNLY","FGFBP2","GZMK","LYAR","TXNIP","GZMM","KLRG1","FCRL6")
breaks <- c(-1.5,0,2)
p <- DotPlot(data,features = heat_genes, dot.scale = 4)+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1), axis.text = element_text(size=12,color="black"),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 12),legend.text = element_text(size = 12))+scale_colour_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")),breaks = breaks, labels = format(breaks))+ scale_size(range = c(4, 8))
ggsave(p, file="PBMC_marker_dotplot_v2.pdf",width=12,height=6)

#remove 9,16
data <- subset(data,idents=c(0:8,10:15,17))
data <- FindNeighbors(data, dims = 1:30, reduction = "harmony", verbose = FALSE)
data <- FindClusters(data, reduction.type = "harmony",resolution = seq(0.1, 1, 0.1),dims.use = 1:30, verbose = FALSE)
data <- RunUMAP(data, dims = 1:30,reduction = "harmony",  verbose = FALSE)
data <- RunTSNE(data, dims = 1:30,reduction = "harmony")
saveRDS(data,file="Tcell_with_coordinate_PC30_v3.rds")

head(data@meta.data)
sapply(grep("^SCT_snn_res",colnames(data@meta.data),value = TRUE),
       function(x) length(unique(data@meta.data[,x])))

data <- SetIdent(data, value="SCT_snn_res.0.9")
p1 <- DimPlot(data, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object = data, reduction = 'umap', label=TRUE, label.size=5)
ggsave(plot_grid(p1,p2),file="data_umap.pdf",width=14,height=6)
p <- DimPlot(object = data, reduction = 'umap', label=TRUE, label.size=5,split.by="orig.ident",ncol=10)
ggsave(p, file="data_umap_split.pdf",width=40,height=5)
p1 <- DimPlot(data, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(object = data, reduction = 'tsne', label=TRUE, label.size=5)
ggsave(plot_grid(p1,p2),file="data_tsne.pdf",width=14,height=6)
p <- DimPlot(object = data, reduction = 'tsne', label=TRUE, label.size=5,split.by="orig.ident",ncol=10)
ggsave(p, file="data_tsne_split.pdf",width=40,height=5)


for (i in levels(data)){
  marker <- FindMarkers(data,ident.1=i, ident.2=NULL, test.use="roc",only.pos=TRUE)
  write.table(marker,file=paste0(i,"_marker.xls"),sep="\t",quote=FALSE)
}


genes <- c("CD3D","CD3E","CD3G","CD8A","CD8B","CD4","CD40LG","CST3","CD68","AIF1","GZMA","KLRD1","NKG7","KLRB1","KLRC1")
breaks <- c(-1,0,1,2)
p <- DotPlot(data,features = genes)+coord_flip()+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=12,color="black"),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 12),legend.text = element_text(size = 12))+scale_colour_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")),breaks = breaks, labels = format(breaks))+ scale_size(range = c(4, 8))
ggsave(p, file="marker_dotplot.pdf",width=13,height=5)


heat_genes <- c("CD4","TCF7","CCR7","SELL","LEF1","CD28","ITGB1","S100A4","LGALS1","KLRB1","GPR183","DNAJA1","DNAJB1","HSP90AB1","HSP90AA1","FOXP3","TIGIT","CTLA4","CD8A","CD8B","PRF1","GZMB","GNLY","FGFBP2","GZMK","LYAR","TXNIP","GZMM","KLRG1","FCRL6")
breaks <- c(-1.5,0,2)
p <- DotPlot(data,features = heat_genes, dot.scale = 4)+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1), axis.text = element_text(size=12,color="black"),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 12),legend.text = element_text(size = 12))+scale_colour_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")),breaks = breaks, labels = format(breaks))+ scale_size(range = c(4, 8))
ggsave(p, file="PBMC_marker_dotplot_v2.pdf",width=12,height=6)


data@meta.data$celltype <- "NA"
data@meta.data[which(data@meta.data$SCT_snn_res.0.9 %in% c(1,3,5,10,11,12,13,14)),]$celltype <- "CD8"
data@meta.data[which(data@meta.data$SCT_snn_res.0.9 %in% c(0,2,4,6,7,8,9,15,16)),]$celltype <- "CD4"
data <- SetIdent(data,value="celltype")

data@meta.data$anno <- "NA"
data@meta.data[which(data@meta.data$SCT_snn_res.0.9 %in% c(0,8,16)),]$anno <- "CD4_Effector_memory"
data@meta.data[which(data@meta.data$SCT_snn_res.0.9 %in% c(2,4)),]$anno <- "CD4_Naive"
data@meta.data[which(data@meta.data$SCT_snn_res.0.9 %in% c(6,7,15)),]$anno <- "CD4_Stress"
data@meta.data[which(data@meta.data$SCT_snn_res.0.9 %in% c(9)),]$anno <- "Treg"
data@meta.data[which(data@meta.data$SCT_snn_res.0.9 %in% c(11)),]$anno <- "CD8_Naive"
data@meta.data[which(data@meta.data$SCT_snn_res.0.9 %in% c(3,12,13)),]$anno <- "CD8_Effector_memory"
data@meta.data[which(data@meta.data$SCT_snn_res.0.9 %in% c(1,5,10,14)),]$anno <- "CD8_Cytotoxicity"
data <- SetIdent(data,value="anno")
saveRDS(data,file="PBMC_Tcell_anno.rds")


for (i in levels(data)){
	marker <- FindMarkers(data,ident.1=i, ident.2=NULL, test.use="roc",only.pos=TRUE)
	write.table(marker,file=paste0(i,"_marker.xls"),sep="\t",quote=FALSE)
}

levels(data) <- c("CD4_Naive","CD4_Effector_memory","CD4_Stress","Treg","CD8_Naive","CD8_Cytotoxicity","CD8_Effector_memory")

p2 <- DimPlot(object = data, reduction = 'tsne')+scale_colour_manual(values = c("#BC3C29FF","#0072B5FF","#B1746FFF","#E18727FF", "#20854EFF","#7876B1FF","#4DBBD5B2"))+xlab("tSNE-1")+ylab("tSNE-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file="PBMC_tsne_anno.pdf",width=5,height=3)

p2 <- DimPlot(object = data, reduction = 'umap')+scale_colour_manual(values = c("#BC3C29FF","#0072B5FF","#B1746FFF","#E18727FF", "#20854EFF","#7876B1FF","#4DBBD5B2"))+xlab("UMAP-1")+ylab("UMAP-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file="PBMC_umap_anno.pdf",width=5,height=3)

heat_genes <- c("CD4","CD40LG","TCF7","CCR7","SELL","LEF1","ITGB1","S100A4","LGALS1","TNFRSF4","DNAJB1","HSP90AA1","HSPA1B","FOXP3","TIGIT","CTLA4","IL2RA","PRF1","GZMB","GNLY","FGFBP2","GZMK","LYAR","TXNIP","GZMM","KLRG1")

levels(data) <- c("CD8_Effector_memory","CD8_Cytotoxicity","CD8_Naive","Treg","CD4_Stress","CD4_Effector_memory","CD4_Naive")

breaks <- c(-1.5,0,2)
p <- DotPlot(data,features = heat_genes, dot.scale = 4)+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=12,color="black"),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 12),legend.text = element_text(size = 12))+scale_colour_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")),breaks = breaks, labels = format(breaks))+ scale_size(range = c(4, 8))
ggsave(p, file="marker_dotplot.pdf",width=12,height=4)



#percent
#Fraction of ecah cell type
num_matrix=matrix(nrow=length(levels(data)),ncol=length(unique(data@meta.data$orig.ident)))
j <- 0
for(i in levels(data)){
    j <- j+1
    cluster=WhichCells(data, idents=i)
    a <- 0
    for (sample in unique(data@meta.data$orig.ident)){
        sample2=rownames(data@meta.data)[which(data@meta.data$orig.ident==sample)]
        num=length(intersect(cluster,sample2))
        a <- a+1
        num_matrix[j,a] <- num
    }
}
colnames(num_matrix)=unique(data@meta.data$orig.ident)
num_matrix2 <- rbind(num_matrix,as.vector(colSums(num_matrix)))
rownames(num_matrix2)=c(levels(data),"Total")
num_matrix2
num_matrix3 <- as.data.frame(num_matrix2)
num_matrix3$total <- rowSums(num_matrix2)
num_matrix3
write.table(num_matrix3, file="Number_of_eachcell_in_eachsample_split.xls",sep="\t",quote=FALSE)

percent_matrix=matrix(nrow=length(levels(data)),ncol=length(unique(data@meta.data$orig.ident)))
j <- 0
for(i in levels(data)){
    j <- j+1
    cluster=WhichCells(data, idents=i)
    a <- 0
    for (sample in unique(data@meta.data$orig.ident)){
        sample2=rownames(data@meta.data)[which(data@meta.data$orig.ident==sample)]
        percent=length(intersect(cluster,sample2))/length(sample2)
        a <- a+1
        percent_matrix[j,a] <- percent
    }
}
colnames(percent_matrix)=unique(data@meta.data$orig.ident)
rownames(percent_matrix)=c(levels(data))
percent_matrix
write.table(percent_matrix, file="Percent_of_eachcell_in_eachsample_split.xls",sep="\t",quote=FALSE)


library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggpubr)

file1="Percent_of_eachcell_in_eachsample_split.xls"
percent1 <- read.csv(file1,header=T,row.names=1,sep="\t")
df <- data.frame(sample = rep(colnames(percent1), each = 7), celltype = rep(rownames(percent1),9), Fraction = as.numeric(as.matrix(percent1)))

pdf("Percent_stackbar_PBMC.pdf",width=8.7,height=3.9)
ggplot(df,aes(fill=factor(celltype,levels=c("CD4_Naive","CD4_Effector_memory","CD4_Stress","Treg","CD8_Naive","CD8_Cytotoxicity","CD8_Effector_memory")),y=Fraction,x=sample))+geom_bar(position="stack", stat="identity")+scale_fill_manual(values=c("#BC3C29FF","#0072B5FF","#B1746FFF","#E18727FF", "#20854EFF","#7876B1FF","#4DBBD5B2"))+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1), axis.text = element_text(size=12,color="black"),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 12),legend.text = element_text(size = 12))
dev.off()

data <- SetIdent(data,value="celltype")
CD4 <- subset(data,ident="CD4")
CD4 <- SetIdent(CD4,value="anno")
CD8 <- subset(data,ident="CD8")
CD8 <- SetIdent(CD8,value="anno")
saveRDS(CD4,file="CD4.anno.rds")
saveRDS(CD8,file="CD8.anno.rds")

levels(CD8) <- c("CD8_Naive","CD8_Cytotoxicity","CD8_Effector_memory")
p2 <- DimPlot(object = CD8, reduction = 'tsne')+scale_colour_manual(values = c("#20854EFF","#7876B1FF","#4DBBD5B2"))+xlab("tSNE-1")+ylab("tSNE-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file="CD8_tsne_anno.pdf",width=5,height=3)

p2 <- DimPlot(object = CD8, reduction = 'umap')+scale_colour_manual(values = c("#20854EFF","#7876B1FF","#4DBBD5B2"))+xlab("UMAP-1")+ylab("UMAP-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file="CD8_umap_anno.pdf",width=5,height=3)

levels(CD4) <- c("CD4_Naive","CD4_Effector_memory","CD4_Stress","Treg")
p2 <- DimPlot(object = CD4, reduction = 'tsne')+scale_colour_manual(values = c("#BC3C29FF","#0072B5FF","#B1746FFF","#E18727FF"))+xlab("tSNE-1")+ylab("tSNE-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file="CD4_tsne_anno.pdf",width=5,height=3)

p2 <- DimPlot(object = CD4, reduction = 'umap')+scale_colour_manual(values = c("#BC3C29FF","#0072B5FF","#B1746FFF","#E18727FF"))+xlab("UMAP-1")+ylab("UMAP-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file="CD4_umap_anno.pdf",width=5,height=3)


