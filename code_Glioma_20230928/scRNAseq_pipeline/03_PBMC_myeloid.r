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


setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/PBMC")

Myeloid <- readRDS("../../RDS_file/Myeloid_cell.rds")

head(Myeloid@meta.data)

Myeloid <- SetIdent(Myeloid,value="group")

PBMC <- subset(Myeloid,ident="Glioma")
PBMC <- SetIdent(PBMC,value="anno")

counts <- PBMC@assays$RNA@counts

data <- CreateSeuratObject(counts = counts, project = "PBM", min.cells = 3, min.features = 200)

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
ribo_genes <- rownames(data@assays$RNA@data)[grep("^RP[SL]",rownames(data@assays$RNA@data))]
data[["percent.ribo"]] <- PercentageFeatureSet(data, features=ribo_genes)
hb_genes <- rownames(data@assays$RNA@data)[grep("^HB[^(P)]",rownames(data@assays$RNA@data))]
data[["percent.hb"]] <- PercentageFeatureSet(data, features=hb_genes)
HSP_genes <- rownames(data)[grep("HSP",rownames(data))]
data[["percent.hsp"]] <- PercentageFeatureSet(data, features=HSP_genes)


data <- SetIdent(data,value="orig.ident")
p <- VlnPlot(object = data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb","percent.hsp"), group.by="orig.ident",  ncol = 5, pt.size=0)
ggsave(p, file="data_sample_QC.pdf",width=40,height=5)


#options(future.globals.maxSize = 8000 * 1024^2)

data<- SCTransform(data, method = "glmGamPoi", assay = "RNA",new.assay.name = "SCT", vars.to.regress = c("nCount_RNA","percent.mt"), return.only.var.genes=FALSE)

data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000, assay="RNA")
saveRDS(data,file="data_sctrans.rds")

DefaultAssay(data) <- "SCT"
length(VariableFeatures(data))

data <- RunPCA(data, features = VariableFeatures(data))
#data <- CellCycleScoring(object = data, s.features = s.genes, g2m.features = g2m.genes, set.ident=TRUE)
p <- PCAPlot(object = data)
ggsave(p,file="data_PCA.pdf",width=5,height=5)

#harmony remove batch effect
p1 <- DimPlot(object = data, reduction = "pca", group.by = "orig.ident")
p2 <- VlnPlot(object = data, features = "PC_1", group.by = "orig.ident")
#plot_grid(p1,p2)
ggsave(plot_grid(p1,p2),file="All_sample_pca_before_harmony.pdf",width=14,height=7)
pdf("RunHarmony.pdf")
data <- RunHarmony(object=data, group.by.vars="orig.ident",assay.use="SCT", plot_convergence = TRUE)
dev.off()
data_harmony_embeddings <- Embeddings(data, 'harmony')
data_harmony_embeddings[1:5, 1:5]
p1 <- DimPlot(object = data, reduction = "harmony", group.by = "orig.ident")
p2 <- VlnPlot(object = data, features = "harmony_1", group.by = "orig.ident")
#plot_grid(p1,p2)
ggsave(plot_grid(p1,p2),file="All_sample_pca_after_harmony_group_by_sample.pdf",width=16,height=8)
p <- ElbowPlot(data,reduction="harmony",ndims=50)
ggsave(p, file="ElbowPlot.pdf",width=5,height=5)

data <- FindNeighbors(data, dims = 1:30, reduction = "harmony", verbose = FALSE)
data <- FindClusters(data, reduction.type = "harmony",resolution = seq(0.1, 1.0, 0.1),dims.use = 1:30, verbose = FALSE)
data <- RunUMAP(data, dims = 1:30,reduction = "harmony",verbose = FALSE)
data <- RunTSNE(data, dims = 1:30,reduction = "harmony",verbose = FALSE)
saveRDS(data,file="PBMC_data.rds")

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


genes <- c("CD3D","CD3E","CD3G","CD8A","CD8B","CD4","CD40LG","CST3","CD68","AIF1","GZMA","KLRD1","NKG7","KLRB1","KLRC1","C1QA","LYZ","VCAN")
p <- DoHeatmap(data,features = genes,size=5,draw.lines=TRUE,combin=TRUE)+scale_fill_gradientn(colors = c("blue", "white", "red"),na.value="white")+theme(text = element_text(size=11),axis.text = element_text(color="black"))
ggsave(p, file="marker_heatmap.pdf",width=20,height=5)

for (i in levels(data)){
  marker <- FindMarkers(data,ident.1=i, ident.2=NULL, test.use="roc",only.pos=TRUE)
  write.table(marker,file=paste0(i,"_marker.xls"),sep="\t",quote=FALSE)
}


genes <- c("LYZ","MNDA","S100A9","S100A8","S100A6","FCN1","VCAN","CD1C","CD1E","CLEC10A","FCER1A","ITGA4","MRC1","TGFBI","THBD","CX3CR1","TMEM119","ITGAX","P2RY12","BIN1","NAV3","CD74","HLA-DPA1","HLA-DPB1","HLA-DRA","RNASE1","MIF","CTSB","CTSL","ENO1","IFI27","IFI6","LY6E","MX1","ISG15","CCL4","CCL4L2","CCL3","CCL3L1","HSPA6","HSPA1A","HSPA1B","DNAJB1")
p <- DoHeatmap(data,features = genes,size=5,draw.lines=TRUE,combin=TRUE)+scale_fill_gradientn(colors = c("blue", "white", "red"),na.value="white")+theme(text = element_text(size=11),axis.text = element_text(color="black"))
ggsave(p, file="marker_heatmap_v2.pdf",width=20,height=5)


#remove 8,16,17,18
#remove doublets and re-cluster
data1 <- subset(data,idents=c(0:7,9:15))
counts <- data1@assays$RNA@counts
data <- CreateSeuratObject(counts = counts, project = "GBM", min.cells = 3, min.features = 200)

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
ribo_genes <- rownames(data@assays$RNA@data)[grep("^RP[SL]",rownames(data@assays$RNA@data))]
data[["percent.ribo"]] <- PercentageFeatureSet(data, features=ribo_genes)
hb_genes <- rownames(data@assays$RNA@data)[grep("^HB[^(P)]",rownames(data@assays$RNA@data))]
data[["percent.hb"]] <- PercentageFeatureSet(data, features=hb_genes)

p1 <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo","percent.hb"), ncol = 5)
ggsave(p1,file="Myeloid_QC.pdf",width=25,height=8)


data <- SCTransform(data, method = "glmGamPoi", assay = "RNA",new.assay.name = "SCT", vars.to.regress = c("nCount_RNA","percent.mt"), return.only.var.genes=FALSE)
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000, assay="RNA")
saveRDS(data,file="Myeloid_sctrans.rds")

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

#30 ok
data <- FindNeighbors(data, dims = 1:30, reduction = "harmony", verbose = FALSE)
data <- FindClusters(data, reduction.type = "harmony",resolution = seq(0.1, 1, 0.1),dims.use = 1:30, verbose = FALSE)
data <- RunUMAP(data, dims = 1:30,reduction = "harmony",  verbose = FALSE)
data <- RunTSNE(data, dims = 1:30,reduction = "harmony")
saveRDS(data,file="Myeloid_with_coordinate_PC30.rds")

head(data@meta.data)
sapply(grep("^SCT_snn_res",colnames(data@meta.data),value = TRUE),
       function(x) length(unique(data@meta.data[,x])))

data <- SetIdent(data, value="SCT_snn_res.0.3")
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


data@meta.data$anno <- "NA"
data@meta.data[which(data@meta.data$SCT_snn_res.0.2 %in% c(0,1,2,4,6)),]$anno <- "Monocyte"
data@meta.data[which(data@meta.data$SCT_snn_res.0.2 %in% c(5)),]$anno <- "DC"
data@meta.data[which(data@meta.data$SCT_snn_res.0.2 %in% c(3)),]$anno <- "MDM"
data@meta.data$celltype <- "NA"
data@meta.data[which(data@meta.data$SCT_snn_res.0.2 %in% c(0,1,2,4,6)),]$celltype <- "Monocyte"
data@meta.data[which(data@meta.data$SCT_snn_res.0.2 %in% c(5)),]$celltype <- "DC"
data@meta.data[which(data@meta.data$SCT_snn_res.0.2 %in% c(3)),]$celltype <- "MDM"

data <- SetIdent(data,value="anno")

saveRDS(data,file="PBMC_myeloid_anno.rds")


levels(data) <- c("Monocyte","DC","MDM")

p2 <- DimPlot(object = data, reduction = 'tsne', pt.size=0.001)+scale_colour_manual(values = c("#bfa2cb","#fbef48",'#0082b8'))+xlab("tSNE-1")+ylab("tSNE-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file="PBMC_tsne_anno.pdf",width=4.2,height=3)

p2 <- DimPlot(object = data, reduction = 'umap', pt.size=0.001)+scale_colour_manual(values = c("#bfa2cb","#fbef48",'#0082b8'))+xlab("UMAP-1")+ylab("UMAP-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file="PBMC_umap_anno.pdf",width=4.2,height=3)

p2 <- DimPlot(object = data, reduction = 'tsne', pt.size=0.001)+scale_colour_manual(values = c("#bfa2cb","#fbef48","#00A0B7B2"))+xlab("tSNE-1")+ylab("tSNE-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file="PBMC_tsne_celltype.pdf",width=4.2,height=3)

p2 <- DimPlot(object = data, reduction = 'umap', pt.size=0.001)+scale_colour_manual(values = c("#bfa2cb","#fbef48","#00A0B7B2"))+xlab("UMAP-1")+ylab("UMAP-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file="PBMC_umap_celltype.pdf",width=4.2,height=3)


#marker
levels(data) <- c("Monocyte","DC","MDM")

genes <- c("LYZ","MNDA","S100A9","S100A8","S100A6","FCN1","VCAN","CD1C","CD1E","CLEC10A","FCER1A","ITGA4","TGFBI","THBD")

Monocyte_markers <- c("LYZ","MNDA","S100A9","S100A8","S100A6","FCN1","VCAN")
DC_markers <- c("CD1C","CD1E","CLEC10A","FCER1A")
MDM_markers <- c("TGFBI","THBD")

features <- list("Monocyte" = Monocyte_markers, "DC" = DC_markers, "MDM" = MDM_markers)
breaks <- c(-2,-1,0,1,2)
p <- DotPlot(data,features = features, dot.scale = 4)+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1), axis.text = element_text(size=12,color="black"),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 12),legend.text = element_text(size = 12))+scale_colour_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")),breaks = breaks, labels = format(breaks))+ scale_size(range = c(4, 8))
ggsave(p, file="marker_dotplot_split_PBMC.pdf",width=8,height=4)


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

#percent of each cell type in all cells
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggpubr)

file1="Percent_of_eachcell_in_eachsample_split.xls"
percent1 <- read.csv(file1,header=T,row.names=1,sep="\t")
df <- data.frame(sample = rep(colnames(percent1), each = 3), celltype = rep(rownames(percent1),9), Fraction = as.numeric(as.matrix(percent1)))

pdf("Percent_stackbar_PBMC.pdf",width=7.9,height=4)
ggplot(df,aes(fill=factor(celltype,levels=c("Monocyte","DC","MDM")),y=Fraction,x=sample))+geom_bar(position="stack", stat="identity")+scale_fill_manual(values=c("#bfa2cb","#fbef48",'#0082b8'))+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1), axis.text = element_text(size=12,color="black"),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 12),legend.text = element_text(size = 12))
dev.off()
