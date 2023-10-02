library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(Matrix)
library(cowplot)
library(dplyr)
library(harmony)
library(reshape2)

setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/cellphonedb/files")

GBM_Tcell <- readRDS("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Tcell/Glioma/GBM_Tcell_anno.rds")

PBMC_Tcell <- readRDS("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Tcell/PBMC/PBMC_Tcell_anno.rds")

Myeloid_GBM <- readRDS("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/Myeloid_anno.rds")

Myeloid_PBMC <- readRDS("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/PBMC/PBMC_myeloid_anno.rds")

All <- readRDS("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/data_anno.rds")
B_vs_NK <- subset(All,idents=c("B","NK"))
B_vs_NK <- SetIdent(B_vs_NK,value="group")
GBM_B_vs_NK <- subset(B_vs_NK,ident="GBM")
PBMC_B_vs_NK <- subset(B_vs_NK,ident="PBMC")


GBM_counts <- All@assays$RNA@data[,c(WhichCells(Myeloid_GBM),WhichCells(GBM_Tcell),WhichCells(GBM_B_vs_NK))]
colnames(GBM_counts) <- gsub("\\.","-",colnames(GBM_counts))
write.table(GBM_counts,file="GBM_cells_count.txt",sep="\t",quote=FALSE)

PBMC_counts <- All@assays$RNA@data[,c(WhichCells(Myeloid_PBMC),WhichCells(PBMC_Tcell),WhichCells(PBMC_B_vs_NK))]
colnames(PBMC_counts) <- gsub("\\.","-",colnames(PBMC_counts))
write.table(PBMC_counts,file="PBMC_cells_count.txt",sep="\t",quote=FALSE)


GBM_meta <- data.frame(Cell=colnames(GBM_counts),celltype=c(Myeloid_GBM@meta.data$anno,GBM_Tcell@meta.data$anno,GBM_B_vs_NK@meta.data$anno))
write.table(GBM_meta,file="GBM_meta.txt",sep="\t",row.names=FALSE,quote=FALSE)


PBMC_meta <- data.frame(Cell=colnames(PBMC_counts),celltype=c(Myeloid_PBMC@meta.data$anno,PBMC_Tcell@meta.data$anno,PBMC_B_vs_NK@meta.data$anno))
write.table(PBMC_meta,file="PBMC_meta.txt",sep="\t",row.names=FALSE,quote=FALSE)


GBM_meta <- data.frame(Cell=colnames(GBM_counts),celltype=c(Myeloid_GBM@meta.data$anno,GBM_Tcell@meta.data$celltype,GBM_B_vs_NK@meta.data$anno))
write.table(GBM_meta,file="GBM_meta_v2.txt",sep="\t",row.names=FALSE,quote=FALSE)

PBMC_meta <- data.frame(Cell=colnames(PBMC_counts),celltype=c(Myeloid_PBMC@meta.data$anno,PBMC_Tcell@meta.data$celltype,PBMC_B_vs_NK@meta.data$anno))
write.table(PBMC_meta,file="PBMC_meta_v2.txt",sep="\t",row.names=FALSE,quote=FALSE)

GBM_meta <- data.frame(Cell=colnames(GBM_counts),celltype=c(Myeloid_GBM@meta.data$celltype,GBM_Tcell@meta.data$celltype,GBM_B_vs_NK@meta.data$anno))
write.table(GBM_meta,file="GBM_meta_v3.txt",sep="\t",row.names=FALSE,quote=FALSE)

PBMC_meta <- data.frame(Cell=colnames(PBMC_counts),celltype=c(Myeloid_PBMC@meta.data$anno,PBMC_Tcell@meta.data$celltype,PBMC_B_vs_NK@meta.data$anno))
write.table(PBMC_meta,file="PBMC_meta_v3.txt",sep="\t",row.names=FALSE,quote=FALSE)


GBM_meta <- data.frame(Cell=colnames(GBM_counts),celltype=c(Myeloid_GBM@meta.data$celltype,GBM_Tcell@meta.data$anno,GBM_B_vs_NK@meta.data$anno))
write.table(GBM_meta,file="GBM_meta_v4.txt",sep="\t",row.names=FALSE,quote=FALSE)

PBMC_meta <- data.frame(Cell=colnames(PBMC_counts),celltype=c(Myeloid_PBMC@meta.data$anno,PBMC_Tcell@meta.data$anno,PBMC_B_vs_NK@meta.data$anno))
write.table(PBMC_meta,file="PBMC_meta_v4.txt",sep="\t",row.names=FALSE,quote=FALSE)
