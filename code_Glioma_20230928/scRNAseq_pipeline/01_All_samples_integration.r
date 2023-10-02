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

dir="/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/result/"
samples <- c(paste0("GBM",c(1:9)),paste0("PBMC",c(1:9)))
sceList <- lapply(samples, function(x){
    #x <- samples[1]
    data <- readRDS(paste0(dir,x,"/",x,"_singlet.rds"))
})

merged <- merge(sceList[[1]], y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],sceList[[6]],sceList[[7]],sceList[[8]],sceList[[9]],sceList[[10]],sceList[[11]],sceList[[12]],sceList[[13]],sceList[[14]],sceList[[15]],sceList[[16]],sceList[[17]],sceList[[18]]), add.cell.ids = samples, project = "merged")


merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
ribo_genes <- rownames(merged@assays$RNA@data)[grep("^RP[SL]",rownames(merged@assays$RNA@data))]
merged[["percent.ribo"]] <- PercentageFeatureSet(merged, features=ribo_genes)
hb_genes <- rownames(merged@assays$RNA@data)[grep("^HB[^(P)]",rownames(merged@assays$RNA@data))]
merged[["percent.hb"]] <- PercentageFeatureSet(merged, features=hb_genes)
HSP_genes <- rownames(merged)[grep("HSP",rownames(merged))]
merged[["percent.hsp"]] <- PercentageFeatureSet(merged, features=HSP_genes)


merged <- SetIdent(merged,value="orig.ident")
p <- VlnPlot(object = merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb","percent.hsp"), group.by="orig.ident",  ncol = 5, pt.size=0)
ggsave(p, file="Merged_sample_QC.pdf",width=40,height=5)


#options(future.globals.maxSize = 8000 * 1024^2)

merged<- SCTransform(merged, method = "glmGamPoi", assay = "RNA",new.assay.name = "SCT", vars.to.regress = c("nCount_RNA","percent.mt"), return.only.var.genes=FALSE)

merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000, assay="RNA")
saveRDS(merged,file="merged_sctrans.rds")

DefaultAssay(merged) <- "SCT"
length(VariableFeatures(merged))
merged <- RunPCA(merged, features = VariableFeatures(merged))
p <- PCAPlot(object = merged)
ggsave(p,file="merged_cellcylce.pdf",width=5,height=5)

#harmony remove batch effect
p1 <- DimPlot(object = merged, reduction = "pca", group.by = "orig.ident")
p2 <- VlnPlot(object = merged, features = "PC_1", group.by = "orig.ident")
#plot_grid(p1,p2)
ggsave(plot_grid(p1,p2),file="All_sample_pca_before_harmony.pdf",width=14,height=7)
pdf("RunHarmony.pdf")
merged <- RunHarmony(object=merged, group.by.vars="orig.ident",assay.use="SCT", plot_convergence = TRUE)
dev.off()
merged_harmony_embeddings <- Embeddings(merged, 'harmony')
merged_harmony_embeddings[1:5, 1:5]
p1 <- DimPlot(object = merged, reduction = "harmony", group.by = "orig.ident")
p2 <- VlnPlot(object = merged, features = "harmony_1", group.by = "orig.ident")
#plot_grid(p1,p2)
ggsave(plot_grid(p1,p2),file="All_sample_pca_after_harmony_group_by_sample.pdf",width=16,height=8)
p <- ElbowPlot(merged,reduction="harmony",ndims=50)
ggsave(p, file="ElbowPlot.pdf",width=5,height=5)

merged <- FindNeighbors(merged, dims = 1:30, reduction = "harmony", verbose = FALSE)
merged <- FindClusters(merged, reduction.type = "harmony",resolution = seq(0.1, 1.0, 0.1),dims.use = 1:30, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:30,reduction = "harmony",verbose = FALSE)
merged <- RunTSNE(merged, dims = 1:30,reduction = "harmony",verbose = FALSE)
saveRDS(merged,file="GBM_merged.rds")


head(merged@meta.data)
sapply(grep("^SCT_snn_res",colnames(merged@meta.data),value = TRUE),
       function(x) length(unique(merged@meta.data[,x])))

merged <- SetIdent(merged, value="SCT_snn_res.0.5")
p1 <- DimPlot(merged, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object = merged, reduction = 'umap', label=TRUE, label.size=5)
ggsave(plot_grid(p1,p2),file="merged_umap.pdf",width=14,height=6)
p <- DimPlot(object = merged, reduction = 'umap', label=TRUE, label.size=5,split.by="orig.ident",ncol=10)
ggsave(p, file="merged_umap_split.pdf",width=40,height=5)
p1 <- DimPlot(merged, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(object = merged, reduction = 'tsne', label=TRUE, label.size=5)
ggsave(plot_grid(p1,p2),file="merged_tsne.pdf",width=14,height=6)
p <- DimPlot(object = merged, reduction = 'tsne', label=TRUE, label.size=5,split.by="orig.ident",ncol=10)
ggsave(p, file="merged_tsne_split.pdf",width=40,height=5)


genes <- c("CD3D","CD3E","CD3G","CD8A","CD8B","CD4","CD40LG","NCAM1","NCR1","KLRF1","CD19","CD79A","CD79B","MS4A1","CD14","LYZ","FCN1", "S100A8","S100A9","MNDA","CD68","ITGAM", "CD163","AIF1", "CST3","C1QA","CD1C","FCER1A","CD1E","CLEC10A","MS4A2","HDC","CTSG","SOX2","SPARC","GFAP","IGFBP7","CLU","FAP","ACTA2","PDGFRA","PDGFRB")
p <- DoHeatmap(merged,features = genes,size=3,draw.lines=TRUE,combin=TRUE)+scale_fill_gradientn(colors = c("blue", "white", "red"),na.value="white")+theme(text = element_text(size=12),axis.text = element_text(color="black"))
ggsave(p, file="marker_heatmap.pdf",width=30,height=12)


for(i in levels(merged)){
    marker <- FindMarkers(merged,ident.1=i,test.use = "roc",only.pos=TRUE,verbose=FALSE)
    write.table(marker,file=paste0("Cluster",i,"_marker.xls"),sep="\t",quote=FALSE)
}

merged <- SetIdent(merged,value="SCT_snn_res.0.5")
merged@meta.data$anno <- "NA"
merged@meta.data$anno[which(merged@meta.data$SCT_snn_res.0.5 %in% c(12,17))] = "NK"
merged@meta.data$anno[which(merged@meta.data$SCT_snn_res.0.5 %in% c(1,2,6,18,21,23,25,28))] = "T"
merged@meta.data$anno[which(merged@meta.data$SCT_snn_res.0.5 %in% c(20))] = "Tumor"
merged@meta.data$anno[which(merged@meta.data$SCT_snn_res.0.5 %in% c(9))] = "B"
merged@meta.data$anno[which(merged@meta.data$SCT_snn_res.0.5 == 13)] = "DC"
merged@meta.data$anno[which(merged@meta.data$SCT_snn_res.0.5 %in% c(3,7,24))] = "Monocyte"
merged@meta.data$anno[which(merged@meta.data$SCT_snn_res.0.5 %in% c(0,4,5,11,14,19,29))] = "Macrophage"
merged@meta.data$anno[which(merged@meta.data$SCT_snn_res.0.5 %in% c(22))] = "Mast"
merged@meta.data$anno[which(merged@meta.data$SCT_snn_res.0.5 %in% c(8,10,15,16,26,27))] = "Doublets"
merged <- SetIdent(merged,value="anno")

merged@meta.data$group <- "NA"
merged@meta.data[grep("GBM",merged@meta.data$orig.ident),]$group <- "GBM"
merged@meta.data[grep("PBMC",merged@meta.data$orig.ident),]$group <- "PBMC"
saveRDS(merged,file="GBM_data_anno.rds")

levels(merged) <- c("Monocyte","Macrophage","DC","Mast","B","NK","T", "Tumor","Doublets")

#remove doublets
data1 <- subset(merged,idents=c("Monocyte","Macrophage","DC","Mast","B","NK","T","Tumor"))
data1 <- SetIdent(data1,value="anno")

data1 <- RunTSNE(data1, dims = 1:30,reduction = "harmony",verbose = FALSE)
data1 <- RunUMAP(data1, dims = 1:30,reduction = "harmony",verbose = FALSE)

head(data1@meta.data)
sapply(grep("^SCT_snn_res",colnames(data1@meta.data),value = TRUE),
       function(x) length(unique(data1@meta.data[,x])))

data1 <- SetIdent(data1, value="SCT_snn_res.0.5")
p1 <- DimPlot(data1, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object = data1, reduction = 'umap', label=TRUE, label.size=5)
ggsave(plot_grid(p1,p2),file="data1_umap.pdf",width=14,height=6)
p <- DimPlot(object = data1, reduction = 'umap', label=TRUE, label.size=5,split.by="orig.ident",ncol=10)
ggsave(p, file="data1_umap_split.pdf",width=40,height=5)
p1 <- DimPlot(data1, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(object = data1, reduction = 'tsne', label=TRUE, label.size=5)
ggsave(plot_grid(p1,p2),file="data1_tsne.pdf",width=14,height=6)
p <- DimPlot(object = data1, reduction = 'tsne', label=TRUE, label.size=5,split.by="orig.ident",ncol=10)
ggsave(p, file="data1_tsne_split.pdf",width=40,height=5)


genes <- c("CD3D","CD3E","CD3G","CD8A","CD8B","CD4","CD40LG","NCAM1","NCR1","KLRF1","CD19","CD79A","CD79B","MS4A1","CD14","LYZ","FCN1", "S100A8","S100A9","MNDA","CD68","ITGAM", "CD163","AIF1", "CST3","C1QA","CD1C","FCER1A","CD1E","CLEC10A","MS4A2","HDC","CTSG","SOX2","SPARC","GFAP","IGFBP7","CLU","FAP","ACTA2","PDGFRA","PDGFRB")
p <- DoHeatmap(data1,features = genes,size=3,draw.lines=TRUE,combin=TRUE)+scale_fill_gradientn(colors = c("blue", "white", "red"),na.value="white")+theme(text = element_text(size=12),axis.text = element_text(color="black"))
ggsave(p, file="marker_heatmap.pdf",width=30,height=12)

#remove tumor
data1 <- SetIdent(data1, value="SCT_snn_res.0.5")
data2 <- subset(data1,idents=c(0:7,9,11:14,17:22,24,25,28,29))
data2 <- SetIdent(data2,value="anno")
saveRDS(data2,file="data_anno.rds")


data <- subset(data2,idents=c("Monocyte","Macrophage","DC","Mast","B","NK","T"))
levels(data) <- c("Monocyte","Macrophage","DC","Mast","B","NK","T")

data@meta.data$orig.ident <- gsub("GBM","Glioma",data@meta.data$orig.ident)
data@meta.data$group <- gsub("GBM","Glioma",data@meta.data$group)
data@meta.data$anno[which(data@meta.data$anno == "Macrophage")] = "Macrophage/Microglia"
data@meta.data$anno[which(data@meta.data$anno == "DC")] = "Dendritic cell"
data <- SetIdent(data,value="anno")

Myeloid <- subset(data,ident=c("Monocyte","Macrophage/Microglia","Dendritic cell"))
saveRDS(Myeloid,file="Myeloid_cell.rds")
infercnv_tumor <- subset(data2,idents=c("T","Tumor"))
saveRDS(infercnv_tumor,file="infercenv_tumor.rds")
T <- subset(data,ident="T")
saveRDS(T,file="T.rds")

data <- SetIdent(data,value="anno")
levels(data) <- c("Monocyte","Macrophage/Microglia","Dendritic cell","Mast","B","NK","T")

saveRDS(data,file="Glioma_vs_PBMC_anno.rds")

p1 <- DimPlot(object = data, reduction = 'tsne', pt.size = 0.001)+scale_colour_manual(values = c("#F39B7FB2","#E64B35B2","#B09C85B2","#7E6148B2","#4DBBD5B2","#00A087B2","#91D1C2B2"))+xlab("tSNE-1")+ylab("tSNE-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.3), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
p2 <- DimPlot(object = data, reduction = 'tsne', pt.size = 0.001, split.by="group")+scale_colour_manual(values = c("#F39B7FB2","#E64B35B2","#B09C85B2","#7E6148B2","#4DBBD5B2","#00A087B2","#91D1C2B2"))+xlab("tSNE-1")+ylab("tSNE-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.3), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p1,file="All_cell_tsne_anno.pdf",width=6,height=4)
ggsave(p2,file="All_cell_tsne_anno_split.pdf",width=8.7,height=4)

p1 <- DimPlot(data, reduction = "tsne", group.by = "orig.ident")+xlab("tSNE-1")+ylab("tSNE-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.3), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p1,file="All_cell_tsne_sample.pdf",width=6.3,height=5)

p2 <- DimPlot(object = data, reduction = 'tsne',group.by="group", pt.size = 0.001)+xlab("tSNE-1")+ylab("tSNE-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.3), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file="All_cell_tsne_group.pdf",width=5.2,height=4)



p2 <- DimPlot(object = data, reduction = 'umap', pt.size = 0.001)+scale_colour_manual(values = c("#F39B7FB2","#E64B35B2","#B09C85B2","#7E6148B2","#4DBBD5B2","#00A087B2","#91D1C2B2"))+xlab("UMAP-1")+ylab("UMAP-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.3), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file="All_cell_umap_anno.pdf",width=6,height=4)

p2 <- DimPlot(object = data, reduction = 'umap', pt.size = 0.001, split.by="group")+scale_colour_manual(values = c("#F39B7FB2","#E64B35B2","#B09C85B2","#7E6148B2","#4DBBD5B2","#00A087B2","#91D1C2B2"))+xlab("UMAP-1")+ylab("UMAP-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.3), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file="All_cell_umap_sample_split.pdf",width=8.7,height=4)

p1 <- DimPlot(data, reduction = "umap", group.by = "orig.ident")+xlab("tSNE-1")+ylab("tSNE-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.3), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p1,file="All_cell_umap_sample.pdf",width=6.3,height=5)

p2 <- DimPlot(object = data, reduction = 'umap',group.by="group", pt.size = 0.001)+xlab("tSNE-1")+ylab("tSNE-2")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.3), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file="All_cell_umap_group.pdf",width=5.2,height=4)

data <- SetIdent(data,value="group")
Glioma <- subset(data,ident="Glioma")
PBMC <- subset(data,ident="PBMC")
p1 <- DimPlot(Glioma, reduction = "umap", group.by = "orig.ident")+xlab("tSNE-1")+ylab("tSNE-2")+scale_colour_manual(values = c("#F8766D","#E88526","#E88526","#B79F00","#93AA00","#5EB300","#00BA38","#00BF74","#00C19F"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.3), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p1,file="Glioma_cell_umap_sample.pdf",width=6.3,height=5)

p2 <- DimPlot(PBMC, reduction = "umap", group.by = "orig.ident")+xlab("tSNE-1")+ylab("tSNE-2")+scale_colour_manual(values = c("#00BFC4","#00B9E3","#00ADFA","#619CFF","#AE87FF","#DB72FB","#F564E3","#FF61C3","#FF699C"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.3), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file="PBMC_cell_umap_sample.pdf",width=6.3,height=5)

p1 <- DimPlot(Glioma, reduction = "tsne", group.by = "orig.ident")+xlab("tSNE-1")+ylab("tSNE-2")+scale_colour_manual(values = c("#F8766D","#E88526","#E88526","#B79F00","#93AA00","#5EB300","#00BA38","#00BF74","#00C19F"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.3), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p1,file="Glioma_cell_tsne_sample.pdf",width=6.3,height=5)

p2 <- DimPlot(PBMC, reduction = "tsne", group.by = "orig.ident")+xlab("tSNE-1")+ylab("tSNE-2")+scale_colour_manual(values = c("#00BFC4","#00B9E3","#00ADFA","#619CFF","#AE87FF","#DB72FB","#F564E3","#FF61C3","#FF699C"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.3), axis.text = element_text(size=10,color="black"),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),legend.text = element_text(size = 10))
ggsave(p2,file="PBMC_cell_tsne_sample.pdf",width=6.3,height=5)

genes <- c("LYZ","FCN1", "S100A8","MNDA","CD68","CD163","C1QA","CD1C","FCER1A","CD1E","CLEC10A","MS4A2","HDC","CTSG","CD19","CD79A","CD79B","MS4A1","NCAM1","NCR1","KLRF1","CD3D","CD3G")

breaks <- c(-2,1,0,1,2)
p <- DotPlot(data,features = genes, dot.scale = 4)+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5), axis.text = element_text(size=12,color="black"),axis.text.x = element_text(size = 12,angle = 45, vjust = 1, hjust=1),axis.text.y = element_text(size = 12),legend.text = element_text(size = 12))+scale_colour_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")),breaks = breaks, labels = format(breaks))+ scale_size(range = c(4, 8))
ggsave(p, file="marker_dotplot_split.pdf",width=10,height=3.5)
ggsave(p, file="marker_dotplot_split_v2.pdf",width=10.5,height=3.3)

data <- SetIdent(data,value="anno")
levels(data) <- c("Monocyte","Macrophage/Microglia","Dendritic cell","Mast","B","NK","T")

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


#
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggpubr)
library(rstatix)

file1="Percent_of_eachcell_in_eachsample_split.xls"
percent1 <- read.csv(file1,header=TRUE,row.names=1,sep="\t")
df <- data.frame(sample = rep(colnames(percent1), each = 7), celltype = rep(rownames(percent1),18), Fraction = as.numeric(as.matrix(percent1)))
df$group <- "NA"
df$group[grep("PBMC",df$sample)] <- "PBMC"
df$group[grep("Glioma",df$sample)] <- "Glioma"
df$group2 <- paste0(df$celltype,"_",df$group)

stat.test <- df %>% group_by(celltype) %>% wilcox_test(Fraction ~ group) %>% adjust_pvalue(method = "bonferroni") %>% add_significance()
stat.test <- stat.test %>% add_xy_position(x = "group")

df$celltype <- factor(df$celltype,levels=c("Monocyte","Macrophage/Microglia","Dendritic cell","Mast","B","NK","T"))
pdf("Percent_bar_PBMC_vs_GBM.pdf",width=10,height=2.5)
bxp <- ggboxplot(df, x = "group", y = "Fraction", fill = "group", palette=c("#D73027","#4575B4"),facet.by = "celltype",nrow=1, scales = "free", add = c("mean_se", "jitter"), order = c("Glioma","PBMC"))+stat_pvalue_manual(stat.test, hide.ns = FALSE, label = "{p.adj.signif}")+scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+theme(title=element_text(size=10,color="black"),axis.text.x=element_text(angle=45,size=10,color="black",hjust=1),axis.text.y=element_text(size=10,color="black"),axis.line=element_line(size=0.3,color="black"),strip.background = element_blank(),plot.title = element_text(size = 10,vjust = 1,hjust = 0.5))+facet_wrap(~celltype, scales = "free", ncol = 7)
bxp
dev.off()
