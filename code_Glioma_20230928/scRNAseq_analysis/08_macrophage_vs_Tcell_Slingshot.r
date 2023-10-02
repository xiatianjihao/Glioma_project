library(Seurat)
library(harmony)
library(slingshot)
library(SingleCellExperiment)
library(reshape2)
library(ggplot2)
library(Matrix)
library(cowplot)
library(dplyr)
library(ggsci)
library(viridis)
library(scales)
library(tradeSeq)
library(grDevices)
library(RColorBrewer)


#MDM
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/trajecotry/MDM")
data <- readRDS("../../../RDS_file/MDM.rds")
data@meta.data$color <- "NA"
data@meta.data$color[which(data@meta.data$anno == "MDM-1")] <- '#0082b8'
data@meta.data$color[which(data@meta.data$anno == "MDM-2")] <- '#3db940'
data@meta.data$color[which(data@meta.data$anno == "MDM-3")] <- '#a4cece'

data <- SetIdent(data,value="anno")
pdf("MDM_umap.pdf",width=5,height=3)
p1 <- DimPlot(data, reduction = "umap")
p1
dev.off()

sce <- as.SingleCellExperiment(data)

sce <- slingshot(sce, clusterLabels = 'anno', start.clus = "MDM-1", reducedDim = 'UMAP')

colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

pdf("slingshot_v1.pdf",width=3.5,height=4)
plot(reducedDims(sce)$UMAP,col = plotcol, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce), lwd=2, col='black')
dev.off()

pdf("slingshotv2.pdf",width=3.5,height=4)
plot(reducedDims(sce)$UMAP,col = plotcol, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce), lwd=2, col='black', type="lineages")
dev.off()


cell_colors <- sce$color
names(cell_colors) <- sce$anno

pdf("slingshot3.pdf",width=3.5,height=4)
plot(reducedDims(sce)$UMAP,col = cell_colors, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black')
dev.off()

pdf("slingshot4.pdf",width=3.5,height=4)
plot(reducedDims(sce)$UMAP,col = cell_colors, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce), lwd = 2, col = 'black')
dev.off()

save(sce,file="sce.rda")

t <- as.data.frame(sce$slingPseudotime_1)
rownames(t) <- colnames(sce)
write.table(t,file="MDM_Slingshot_pseudotime.xls",sep="\t",quote=FALSE)

color.bar <- function(lut, min, max=-min, nticks=5, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
    pdf("MDM_colorbar_pesudotime.pdf", width=2,height=3)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
    dev.off()
}
color.bar(rev(colorRampPalette(brewer.pal(11,'Spectral')[-4])(100)),0,round(summary(sce$slingPseudotime_1)[6]))

#MG
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/trajecotry/MG")

data <- readRDS("../../../RDS_file/MG.rds")
data@meta.data$color <- "NA"
data@meta.data$color[which(data@meta.data$anno == "MG-1")] <- "#9f6164"
data@meta.data$color[which(data@meta.data$anno == "MG-2")] <- "#cc9f66"
data@meta.data$color[which(data@meta.data$anno == "MG-3")] <- "#fbac00"

data <- SetIdent(data,value="anno")
pdf("MG_umap.pdf",width=5,height=3)
p1 <- DimPlot(data, reduction = "umap")
p1
dev.off()

sce <- as.SingleCellExperiment(data)

sce <- slingshot(sce, clusterLabels = 'anno', start.clus = "MG-3", end.clus = c("MG-2","MG-1"), reducedDim = 'UMAP')
getLineages(sce)

summary(sce$slingPseudotime_1)
colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

pdf("slingshot_MG1_v1.pdf",width=3.5,height=4)
plot(reducedDims(sce)$UMAP,col = plotcol, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce), lwd=2, col='black')
dev.off()

pdf("slingshot_MG1_v2.pdf",width=3.5,height=4)
plot(reducedDims(sce)$UMAP,col = plotcol, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce), lwd=2, col='black', type="lineages")
dev.off()

summary(sce$slingPseudotime_2)
colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))
plotcol <- colors[cut(sce$slingPseudotime_2, breaks=100)]

pdf("slingshot_MG2_v1.pdf",width=3.5,height=4)
plot(reducedDims(sce)$UMAP,col = plotcol, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce), lwd=2, col='black')
dev.off()

pdf("slingshot2_MG2_v2.pdf",width=3.5,height=4)
plot(reducedDims(sce)$UMAP,col = plotcol, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce), lwd=2, col='black', type="lineages")
dev.off()

cell_colors <- sce$color
names(cell_colors) <- sce$anno

pdf("slingshot3.pdf",width=3.5,height=4)
plot(reducedDims(sce)$UMAP,col = cell_colors, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black')
dev.off()

pdf("slingshot4.pdf",width=3.5,height=4)
plot(reducedDims(sce)$UMAP,col = cell_colors, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce), lwd = 2, col = 'black')
dev.off()
save(sce,file="sce.rda")

t <- as.data.frame(sce$slingPseudotime_1)
rownames(t) <- colnames(sce)
write.table(t,file="MG_Slingshot_pseudotime1.xls",sep="\t",quote=FALSE)

t <- as.data.frame(sce$slingPseudotime_2)
rownames(t) <- colnames(sce)
write.table(t,file="MG_Slingshot_pseudotime2.xls",sep="\t",quote=FALSE)


color.bar <- function(lut, min, max=-min, nticks=5, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
    pdf("colorbar_pesudotime_MG1.pdf", width=2,height=3)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
    dev.off()
}
color.bar(rev(colorRampPalette(brewer.pal(11,'Spectral')[-4])(100)),0,round(summary(sce$slingPseudotime_1)[6]))


color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
    pdf("colorbar_pesudotime_MG2.pdf", width=2,height=3)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
    dev.off()
}
color.bar(rev(colorRampPalette(brewer.pal(11,'Spectral')[-4])(100)),0,round(summary(sce$slingPseudotime_2)[6]))

#CD8
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Tcell/Glioma/trajectory/CD8")
data <- readRDS("../../CD8.anno.rds")
data@meta.data$color <- "NA"
data@meta.data$color[which(data@meta.data$anno == "CD8_Cytotoxicity")] <- "#7876B1FF"
data@meta.data$color[which(data@meta.data$anno == "CD8_Exhaustion")] <- "#6F99ADFF"
data@meta.data$color[which(data@meta.data$anno == "CD8_Stress")] <- "#725663FF"

data <- SetIdent(data,value="anno")
pdf("CD8_umap.pdf",width=5,height=3)
p1 <- DimPlot(data, reduction = "umap")
p1
dev.off()

sce <- as.SingleCellExperiment(data)
sce <- slingshot(sce, clusterLabels = 'anno',start.clus = "CD8_Cytotoxicity",end.clus = "CD8_Stress",reducedDim = "UMAP",stretch = 3)

summary(sce$slingPseudotime_1)

sce_ex <- sce[, sce$anno %in% c("CD8_Cytotoxicity","CD8_Exhaustion")]
sce_stress <- sce[, sce$anno %in% c("CD8_Cytotoxicity","CD8_Stress")]

colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))
plotcol <- colors[cut(sce_stress$slingPseudotime_1, breaks=100)]

pdf("slingshot_stress_v1.pdf",width=3.5,height=4)
plot(reducedDims(sce_stress)$UMAP,col = plotcol, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce_stress), lwd=2, col='black')
dev.off()

pdf("slingshot_stress_v2.pdf",width=3.5,height=4)
plot(reducedDims(sce_stress)$UMAP,col = plotcol, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce_stress), lwd=2, col='black', type="lineages")
dev.off()

summary(sce$slingPseudotime_2)
colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))
plotcol <- colors[cut(sce_ex$slingPseudotime_2, breaks=100)]


pdf("slingshot_Exhaustion_v1.pdf",width=3.5,height=4)
plot(reducedDims(sce_ex)$UMAP,col = plotcol, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce_ex), lwd=2, col='black')
dev.off()

pdf("slingshot2_Exhaustion_v2.pdf",width=3.5,height=4)
plot(reducedDims(sce_ex)$UMAP,col = plotcol, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce_ex), lwd=2, col='black', type="lineages")
dev.off()

cell_colors <- sce$color
names(cell_colors) <- sce$anno

pdf("slingshot3.pdf",width=3.5,height=4)
plot(reducedDims(sce)$UMAP,col = cell_colors, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black')
dev.off()

pdf("slingshot4.pdf",width=3.5,height=4)
plot(reducedDims(sce)$UMAP,col = cell_colors, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce), lwd = 2, col = 'black')
dev.off()

save(sce,file="sce.rda")

t <- as.data.frame(sce$slingPseudotime_1)
rownames(t) <- colnames(sce)
write.table(t,file="CD8_Slingshot_pseudotime1.xls",sep="\t",quote=FALSE)


t <- as.data.frame(sce$slingPseudotime_2)
rownames(t) <- colnames(sce)
write.table(t,file="CD8_Slingshot_pseudotime2.xls",sep="\t",quote=FALSE)

color.bar <- function(lut, min, max=-min, nticks=5, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
    pdf("colorbar_pesudotime_CD8_Stress.pdf", width=2,height=3)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
    dev.off()
}
color.bar(rev(colorRampPalette(brewer.pal(11,'Spectral')[-4])(100)),0,round(summary(sce$slingPseudotime_1)[6]))

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
    pdf("colorbar_pesudotime_CD8_Exhaution.pdf", width=2,height=3)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
    dev.off()
}
color.bar(rev(colorRampPalette(brewer.pal(11,'Spectral')[-4])(100)),0,round(summary(sce$slingPseudotime_2)[6]))

#CD4
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Tcell/Glioma/trajectory/CD4")
data1 <- readRDS("../../CD4.anno.rds")
data <- subset(data1,idents=c("CD4_Naive","CD4_Memory","CD4_Stress"))
data@meta.data$color <- "NA"
data@meta.data$color[which(data@meta.data$anno == "CD4_Naive")] <- "#BC3C29FF"
data@meta.data$color[which(data@meta.data$anno == "CD4_Memory")] <- "#ADB17DFF"
data@meta.data$color[which(data@meta.data$anno == "CD4_Stress")] <- "#B1746FFF"

data <- SetIdent(data,value="anno")
pdf("CD4_umap.pdf",width=5,height=3)
p1 <- DimPlot(data, reduction = "umap")
p1
dev.off()

sce <- as.SingleCellExperiment(data)
sce <- slingshot(sce, clusterLabels = 'anno', start.clus = "CD4_Naive", reducedDim = 'UMAP')

summary(sce$slingPseudotime_1)

colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-4])(200))
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=200)]

pdf("slingshot.pdf",width=3.5,height=4)
plot(reducedDims(sce)$UMAP,col = plotcol, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce), lwd=2, col='black')
dev.off()

pdf("slingshot2.pdf",width=3.5,height=4)
plot(reducedDims(sce)$UMAP,col = plotcol, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce), lwd=2, col='black', type="lineages")
dev.off()

cell_colors <- sce$color
names(cell_colors) <- sce$anno

pdf("slingshot3.pdf",width=3.5,height=4)
plot(reducedDims(sce)$UMAP,col = cell_colors, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black')
dev.off()

pdf("slingshot4.pdf",width=3.5,height=4)
plot(reducedDims(sce)$UMAP,col = cell_colors, pch=16, asp = 1, cex=0.5)
lines(SlingshotDataSet(sce), lwd = 2, col = 'black')
dev.off()

save(sce,file="sce.rda")

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
    pdf("colorbar_pesudotime_Exhaution.pdf", width=2,height=3)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
    dev.off()
}
color.bar(rev(colorRampPalette(brewer.pal(11,'Spectral')[-4])(100)),0,round(summary(sce$slingPseudotime_1)[6]))

t <- as.data.frame(sce$slingPseudotime_1)
rownames(t) <- colnames(sce)
write.table(t,file="CD4_Slingshot_pseudotime.xls",sep="\t",quote=FALSE)