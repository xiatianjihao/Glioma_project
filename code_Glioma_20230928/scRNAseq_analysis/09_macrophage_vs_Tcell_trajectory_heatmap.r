suppressPackageStartupMessages({
library("reshape2")
library("stringr")
library("plyr")
library("dplyr")
library("ggpubr")
library("ggsci")
library("ggrastr")
library("ggvenn")
library("data.table")
library("sscClust")
library("gam")
library("colorRamps")
library("clusterProfiler")
library("org.Hs.eg.db")
library("SingleCellExperiment")
library("Seurat")
})

doParallel::registerDoParallel(cores=128)
source("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/code_Glioma_project_upload/analysis/source_dyn_func.r")
source("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/code_Glioma_project_upload/analysis/function.r")

#MDM
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/trajecotry/MDM")
stype = "MDM"
#
tfs = read.table("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/code_Glioma_project_upload/analysis/TF.xls", sep="\t",header=F)$V1
tfs = unique(as.vector(tfs))

MDM <- readRDS("../../RDS_file/MDM.rds")
MDM$barcode <- colnames(MDM)
MDM$UMAP_1 <- MDM@reductions$umap@cell.embeddings[,1]
MDM$UMAP_2 <- MDM@reductions$umap@cell.embeddings[,2]

diffmap <- MDM@meta.data
rownames(diffmap) <- diffmap$barcode
diffmap <- diffmap[,c("UMAP_1","UMAP_2")]

# dpt_pseudotime
dpt <- read.csv("MDM_Slingshot_pseudotime.xls",sep="\t")
diffmap$dpt_pseudotime <- dpt[colnames(MDM),1]

DefaultAssay(MDM) <- "RNA"
MDM <- ScaleData(MDM, features = rownames(MDM))
sce <- as.SingleCellExperiment(MDM)
sce <- sce[,rownames(diffmap)]
reducedDim(sce, "diffmap") <- diffmap[,c("UMAP_1","UMAP_2")]
sce$dpt_pseudotime <- diffmap$dpt_pseudotime
sce$meta.cluster <- as.character(sce$anno)

sce <- convertDPTperc(sce)
sce$cluster.name <- sce$anno

colSet <- list()
color <- c('#0082b8','#3db940','#a4cece')
names(color) <- c("MDM-1","MDM-2","MDM-3")
colSet[[1]] <- color

### density
dat <- colData(sce)[,c("dpt_pseudotime","dpt_order_perc","meta.cluster")] %>% as.data.frame
pdf("Density.pdf",width=5,height=3)
p <- ggplot(dat, aes(x=dpt_order_perc, fill=meta.cluster)) +
      geom_density(size=0.2, adjust=1.5, n=50, trim=F, color="black", alpha=0.75) +
      theme_classic2() +
      scale_fill_manual(values=colSet[[1]][unique(dat$meta.cluster)])
print(p)
dens_p <- p + theme_void() + theme(plot.margin=margin(t=0,r=0,b=0,l=0,unit="cm"),legend.position="none")+xlim(0,1)+scale_y_continuous(expand = c(0,0))+scale_x_continuous(expand = c(0,0))
dev.off()

### DEG
assay(sce, "exprs") <- MDM@assays$RNA@scale.data[rownames(assay(sce,"logcounts")),colnames(assay(sce,"logcounts"))]
gam.test <- testProcess(sce, tfs)
saveRDS(gam.test, file="MDM.dyn_genes.rds")

genes <- intersect(tfs, rownames(assay(sce,"exprs")))
a <- assay(sce, "exprs")[genes,]

### heatmap
p.cutoff = 0.05
c.cutoff = 0.2
#
gam.test$sig = ifelse(gam.test$adj.p < p.cutoff & abs(gam.test$coef) > c.cutoff, T, F)
#
plot.gene <- as.character(gam.test[gam.test$sig, "geneSymbol"])
highlight.genes <- c("IFITM3","TNF","IFI16","IL1A","CCL5","CD74","HLA-DPB1","HLA-DRA","HLA-DMA","HLA-DPA1","CD86","IL1A","CCL5","VEGFA","CTSA","MMP14","MIF","ENO1","ALDOA","LDHA","EPAS1")
annot=rowAnnotation(gene=anno_mark(at=match(highlight.genes, plot.gene), labels=highlight.genes))
colSet[['dpt_order_perc']]=colorRamp2(c(0,max(sce$dpt_order_perc)),c("white","darkblue"))
#
set.seed(1)
sce$meta.cluster = as.character(sce$meta.cluster)
sce$cluster.name = as.character(sce$cluster.name)
sce$cluster.name = factor(sce$cluster.name)
pdf("MDM_All_gene.pdf",width=8,height=8)
hm.obj = heatmap_sm(sce[plot.gene,], assay.name="exprs",
                 colSet=colSet,
                 out.prefix="MDM.dyn_genes.heatmap",
                 columns=c("dpt_order_perc"), 
                 ncell.downsample=500,
                 use_raster=T, raster_quality=2,
                 columns.order="dpt_order_perc",
                 show_column_names=F, show_row_names=F,
                 row_names_side="right",
                 row_gap = unit(1, "mm"), column_gap = unit(0, "mm"),
                 border=F, ann.bar.height=0.5, palette.name="blue2green2red",
                 pdf.width=10.5, pdf.height=7, do.scale=F,
                 z.lo=-0.5, z.hi=0.5, z.step=0.5, 
                 do.clustering.row=T, do.clustering.col=F,
                 dend.row=T,  right_annotation=annot,
                 clustering.distance="cosine",
                 clustering.method="ward.D2",
                 #row_split=5
                 row_km=3, row_km_repeats=200, #row_title=NULL
                 )
print(hm.obj)
dev.off()

pdf("MDM_genes.heatmap_addDensity_All_gene_C3.pdf", width=5.5, height=4)
HeatmapAnnotation(ggplot=anno_empty(height=unit(1.2, "cm")), which="column") %v% hm.obj

decorate_annotation("ggplot", {
  vp = current.viewport()$name
  print(dens_p, vp=vp)
})
dev.off()

#Microglia MG1
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Myeloid/Glioma/trajecotry/MG")
seurat <- readRDS("../../RDS_file/MG.rds")
MG <- subset(seurat,idents=c("MG-1","MG-3"))
DefaultAssay(MG) <- "RNA"
MG <- ScaleData(MG, features = rownames(MG))
sce = as.SingleCellExperiment(MG)

MG$barcode <- colnames(MG)
MG$UMAP_1 <- MG@reductions$umap@cell.embeddings[,1]
MG$UMAP_2 <- MG@reductions$umap@cell.embeddings[,2]

diffmap <- MG@meta.data
rownames(diffmap) <- diffmap$barcode
diffmap <- diffmap[,c("UMAP_1","UMAP_2")]
  
# dpt_pseudotime
dpt = read.csv("MG_Slingshot_pseudotime2.xls",sep="\t")
diffmap$dpt_pseudotime = dpt[colnames(MG),1]
 
sce = sce[,rownames(diffmap)]
reducedDim(sce, "diffmap") = diffmap[,c("UMAP_1","UMAP_2")]
sce$dpt_pseudotime = diffmap$dpt_pseudotime
sce$meta.cluster = as.character(sce$anno)
#
sce = convertDPTperc(sce)
sce$cluster.name = sce$anno

####
colSet <- list()
color <- c("#9f6164","#fbac00")
names(color) <- c("MG-1","MG-3")
colSet[[1]] <- color

### density
dat <- colData(sce)[,c("dpt_pseudotime","dpt_order_perc","meta.cluster")] %>% as.data.frame
pdf("Density_MG1.pdf",width=5,height=3)
p <- ggplot(dat, aes(x=dpt_order_perc, fill=meta.cluster)) +
      geom_density(size=0.2, adjust=1.5, n=50, trim=F, color="black", alpha=0.75) +
      theme_classic2() +
      scale_fill_manual(values=colSet[[1]][unique(dat$meta.cluster)])
print(p)
dens_p = p + theme_void() + theme(plot.margin=margin(t=0,r=0,b=0,l=0,unit="cm"),legend.position="none")+xlim(0,1)+scale_y_continuous(expand = c(0,0))+scale_x_continuous(expand = c(0,0))
dev.off()

### DEG
assay(sce, "exprs") <- MG@assays$RNA@scale.data[rownames(assay(sce,"logcounts")),colnames(assay(sce,"logcounts"))]
gam.test <- testProcess(sce, tfs)
saveRDS(gam.test, file="MG.dyn_genes_MG1.rds")

genes <- intersect(tfs, rownames(assay(sce,"exprs")))
a <- assay(sce, "exprs")[genes,]

### heatmap
p.cutoff = 0.05
c.cutoff = 0.2
#
gam.test$sig = ifelse(gam.test$adj.p < p.cutoff & abs(gam.test$coef) > c.cutoff, T, F)
#
plot.gene = as.character(gam.test[gam.test$sig, "geneSymbol"])
highlight.genes = c("IFI27","IFI6","LY6E","IFITM3","ISG15","MX1","CCL4","CCL5","CCL3","CCL3L1","CCL4L2","TMEM119","P2RY12", "TREM2","C1QA","C1QB","C1QC","SPP1","APOE","CX3CR1")
annot=rowAnnotation(gene=anno_mark(at=match(highlight.genes, plot.gene), labels=highlight.genes))
colSet[['dpt_order_perc']]=colorRamp2(c(0,max(sce$dpt_order_perc)),c("white","darkblue"))
#colSet[['cluster.name']] = structure(nam.conv$mcolor, names=sce$cluster.name)
#
set.seed(1)
sce$meta.cluster = as.character(sce$meta.cluster)
sce$cluster.name = as.character(sce$cluster.name)
sce$cluster.name = factor(sce$cluster.name)
pdf("MG_All_gene_v1.pdf",width=8,height=8)
hm.obj = heatmap_sm(sce[plot.gene,], assay.name="exprs",
                 colSet=colSet,
                 out.prefix="MG.dyn_genes.heatmap",
                 columns=c("dpt_order_perc"), 
                 ncell.downsample=500,
                 use_raster=T, raster_quality=2,
                 columns.order="dpt_order_perc",
                 show_column_names=F, show_row_names=F,
                 row_names_side="right",
                 row_gap = unit(1, "mm"), column_gap = unit(0, "mm"),
                 border=F, ann.bar.height=0.5, palette.name="blue2green2red",
                 pdf.width=10.5, pdf.height=7, do.scale=F,
                 z.lo=-0.5, z.hi=0.5, z.step=0.5, 
                 do.clustering.row=T, do.clustering.col=F,
                 dend.row=T,  right_annotation=annot,
                 clustering.distance="cosine",
                 clustering.method="ward.D2",
                 #row_split=5
                 row_km=3, row_km_repeats=200, #row_title=NULL
                 )
print(hm.obj)
dev.off()

pdf("MG_genes.heatmap_addDensity_All_gene_MG1_C3.pdf", width=5.5, height=4)
HeatmapAnnotation(ggplot=anno_empty(height=unit(1.2, "cm")), which="column") %v% hm.obj

decorate_annotation("ggplot", {
  vp = current.viewport()$name
  print(dens_p, vp=vp)
})
dev.off()

#Microglia MG2
seurat <- readRDS("../../RDS_file/MG.rds")
MG <- subset(seurat,idents=c("MG-2","MG-3"))
DefaultAssay(MG) <- "RNA"
MG <- ScaleData(MG, features = rownames(MG))
sce <- as.SingleCellExperiment(MG)

diffmap <- MG@meta.data
rownames(diffmap) <- diffmap$barcode
diffmap <- diffmap[,c("UMAP_1","UMAP_2")]

#dpt_pseudotime
dpt = read.csv("MG_Slingshot_pseudotime2.xls",sep="\t")
diffmap$dpt_pseudotime = dpt[colnames(MG),1]
  
## modified sce
sce = sce[,rownames(diffmap)]
reducedDim(sce, "diffmap") = diffmap[,c("UMAP_1","UMAP_2")]
sce$dpt_pseudotime = diffmap$dpt_pseudotime
sce$meta.cluster = as.character(sce$anno)
#
sce = convertDPTperc(sce)
sce$cluster.name = sce$anno

####
colSet = list()
color = c("#cc9f66","#fbac00")
names(color) <- c("MG-2","MG-3")
colSet[[1]] <- color

### density
dat = colData(sce)[,c("dpt_pseudotime","dpt_order_perc","meta.cluster")] %>% as.data.frame
pdf("Density_MG2.pdf",width=5,height=3)
p = ggplot(dat, aes(x=dpt_order_perc, fill=meta.cluster)) +
      geom_density(size=0.2, adjust=1.5, n=50, trim=F, color="black", alpha=0.75) +
      theme_classic2() +
      scale_fill_manual(values=colSet[[1]][unique(dat$meta.cluster)])
print(p)
dens_p = p + theme_void() + theme(plot.margin=margin(t=0,r=0,b=0,l=0,unit="cm"),legend.position="none")+xlim(0,1)+scale_y_continuous(expand = c(0,0))+scale_x_continuous(expand = c(0,0))
dev.off()

### DEG
assay(sce, "exprs") <- MG@assays$RNA@scale.data[rownames(assay(sce,"logcounts")),colnames(assay(sce,"logcounts"))]
gam.test = testProcess(sce, tfs)
saveRDS(gam.test, file="MG.dyn_genes_MG2.rds")

genes <- intersect(tfs, rownames(assay(sce,"exprs")))
a <- assay(sce, "exprs")[genes,]

### heatmap
p.cutoff = 0.05
c.cutoff = 0.2
#
gam.test$sig = ifelse(gam.test$adj.p < p.cutoff & abs(gam.test$coef) > c.cutoff, T, F)
#
plot.gene = as.character(gam.test[gam.test$sig, "geneSymbol"])
highlight.genes = c("IFI27","IFI6","LY6E","IFITM3","ISG15","MX1","HSPA1A","HSPA1B","DNAJB1","HSPA6","HSP90AA1","HSPH1")
annot=rowAnnotation(gene=anno_mark(at=match(highlight.genes, plot.gene), labels=highlight.genes))
colSet[['dpt_order_perc']]=colorRamp2(c(0,max(sce$dpt_order_perc)),c("white","darkblue"))
#colSet[['cluster.name']] = structure(nam.conv$mcolor, names=sce$cluster.name)
#
set.seed(1)
sce$meta.cluster = as.character(sce$meta.cluster)
sce$cluster.name = as.character(sce$cluster.name)
sce$cluster.name = factor(sce$cluster.name)
pdf("MG_All_gene_MG2.pdf",width=8,height=8)
hm.obj = heatmap_sm(sce[plot.gene,], assay.name="exprs",
                 colSet=colSet,
                 out.prefix="MG.dyn_genes.heatmap",
                 columns=c("dpt_order_perc"), 
                 ncell.downsample=500,
                 use_raster=T, raster_quality=2,
                 columns.order="dpt_order_perc",
                 show_column_names=F, show_row_names=F,
                 row_names_side="right",
                 row_gap = unit(1, "mm"), column_gap = unit(0, "mm"),
                 border=F, ann.bar.height=0.5, palette.name="blue2green2red",
                 pdf.width=10.5, pdf.height=7, do.scale=F,
                 z.lo=-0.5, z.hi=0.5, z.step=0.5, 
                 do.clustering.row=T, do.clustering.col=F,
                 dend.row=T,  right_annotation=annot,
                 clustering.distance="cosine",
                 clustering.method="ward.D2",
                 #row_split=5
                 row_km=2, row_km_repeats=200, #row_title=NULL
                 )
print(hm.obj)
dev.off()

pdf("MG_genes.heatmap_addDensity_All_gene_MG2_C2.pdf", width=5.5, height=4)
HeatmapAnnotation(ggplot=anno_empty(height=unit(1.2, "cm")), which="column") %v% hm.obj

decorate_annotation("ggplot", {
  vp = current.viewport()$name
  print(dens_p, vp=vp)
})
dev.off()

#CD8 Exhaustion
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Tcell/Glioma/trajectory/CD8")

seurat <- readRDS("../../GBM_Tcell_anno.rds")
CD8 <- subset(seurat,idents=c("CD8_Cytotoxicity","CD8_Exhaustion"))
DefaultAssay(CD8) <- "RNA"
CD8 <- ScaleData(CD8, features = rownames(CD8))
sce = as.SingleCellExperiment(CD8)

diffmap <- CD8@meta.data
rownames(diffmap) <- diffmap$barcode
diffmap <- diffmap[,c("UMAP_1","UMAP_2")]

# dpt_pseudotime
dpt = read.csv("CD8_Slingshot_pseudotime2.xls",sep="\t")
diffmap$dpt_pseudotime = dpt[colnames(CD8),1]
  
## modified sce
sce = sce[,rownames(diffmap)]
reducedDim(sce, "diffmap") = diffmap[,c("UMAP_1","UMAP_2")]
sce$dpt_pseudotime = diffmap$dpt_pseudotime
sce$meta.cluster = as.character(sce$anno)
#
sce = convertDPTperc(sce)
sce$cluster.name = sce$anno

####
colSet = list()
color = c("#7876B1FF","#6F99ADFF")
names(color) <- c("CD8_Cytotoxicity","CD8_Exhaustion")
colSet[[1]] <- color

### density
dat = colData(sce)[,c("dpt_pseudotime","dpt_order_perc","meta.cluster")] %>% as.data.frame
pdf("Density_Exhaustion.pdf",width=5,height=3)
p = ggplot(dat, aes(x=dpt_order_perc, fill=meta.cluster)) +
      geom_density(size=0.2, adjust=1.5, n=50, trim=F, color="black", alpha=0.75) +
      theme_classic2() +
      scale_fill_manual(values=colSet[[1]][unique(dat$meta.cluster)])
print(p)
dens_p = p + theme_void() + theme(plot.margin=margin(t=0,r=0,b=0,l=0,unit="cm"),legend.position="none")+xlim(0,1)+scale_y_continuous(expand = c(0,0))+scale_x_continuous(expand = c(0,0))
dev.off()


### DEG
assay(sce, "exprs") <- CD8@assays$RNA@scale.data[rownames(assay(sce,"logcounts")),colnames(assay(sce,"logcounts"))]
gam.test = testProcess(sce, tfs)
saveRDS(gam.test, file="CD8.dyn_genes_Ex.rds")

genes <- intersect(tfs, rownames(assay(sce,"exprs")))
a <- assay(sce, "exprs")[genes,]

### heatmap
p.cutoff = 0.05
c.cutoff = 0.2
#
gam.test$sig = ifelse(gam.test$adj.p < p.cutoff & abs(gam.test$coef) > c.cutoff, T, F)
#
plot.gene = as.character(gam.test[gam.test$sig, "geneSymbol"])
highlight.genes = c("GZMB","PRF1","PDCD1","TIGIT","AHR","NR4A1","TOX","TOX2","ENTPD1","CTLA4","HAVCR2","CXCL13","IFNG","LAG3")
annot=rowAnnotation(gene=anno_mark(at=match(highlight.genes, plot.gene), labels=highlight.genes))
colSet[['dpt_order_perc']]=colorRamp2(c(0,max(sce$dpt_order_perc)),c("white","darkblue"))
#colSet[['cluster.name']] = structure(nam.conv$mcolor, names=sce$cluster.name)
#
set.seed(1)
sce$meta.cluster = as.character(sce$meta.cluster)
sce$cluster.name = as.character(sce$cluster.name)
sce$cluster.name = factor(sce$cluster.name)
pdf("CD8_All_gene_v1.pdf",width=8,height=8)
hm.obj = heatmap_sm(sce[plot.gene,], assay.name="exprs",
                 colSet=colSet,
                 out.prefix="CD8.dyn_genes.heatmap",
                 columns=c("dpt_order_perc"), 
                 ncell.downsample=500,
                 use_raster=T, raster_quality=2,
                 columns.order="dpt_order_perc",
                 show_column_names=F, show_row_names=F,
                 row_names_side="right",
                 row_gap = unit(1, "mm"), column_gap = unit(0, "mm"),
                 border=F, ann.bar.height=0.5, palette.name="blue2green2red",
                 pdf.width=10.5, pdf.height=7, do.scale=F,
                 z.lo=-0.5, z.hi=0.5, z.step=0.5, 
                 do.clustering.row=T, do.clustering.col=F,
                 dend.row=T,  right_annotation=annot,
                 clustering.distance="cosine",
                 clustering.method="ward.D2",
                 #row_split=5
                 row_km=4, row_km_repeats=200, #row_title=NULL
                 )
print(hm.obj)
dev.off()


pdf("CD8_genes.heatmap_addDensity_All_gene_Exhaustion.pdf", width=5.5, height=4)
HeatmapAnnotation(ggplot=anno_empty(height=unit(1.2, "cm")), which="column") %v% hm.obj

decorate_annotation("ggplot", {
  vp = current.viewport()$name
  print(dens_p, vp=vp)
})
dev.off()

##CD8 Stress
seurat <- readRDS("../../GBM_Tcell_anno.rds")
CD8 <- subset(seurat,idents=c("CD8_Cytotoxicity","CD8_Stress"))
DefaultAssay(CD8) <- "RNA"
CD8 <- ScaleData(CD8, features = rownames(CD8))
sce = as.SingleCellExperiment(CD8)

diffmap <- CD8@meta.data
rownames(diffmap) <- diffmap$barcode
diffmap <- diffmap[,c("UMAP_1","UMAP_2")]

# dpt_pseudotime
dpt = read.csv("CD8_Slingshot_pseudotime1.xls",sep="\t")
diffmap$dpt_pseudotime = dpt[colnames(CD8),1]
  
## modified sce
sce = sce[,rownames(diffmap)]
reducedDim(sce, "diffmap") = diffmap[,c("UMAP_1","UMAP_2")]
sce$dpt_pseudotime = diffmap$dpt_pseudotime
sce$meta.cluster = as.character(sce$anno)
#
sce = convertDPTperc(sce)
sce$cluster.name = sce$anno

####
colSet = list()
color = c("#7876B1FF","#725663FF")
names(color) <- c("CD8_Cytotoxicity","CD8_Stress")
colSet[[1]] <- color

### density
dat = colData(sce)[,c("dpt_pseudotime","dpt_order_perc","meta.cluster")] %>% as.data.frame
pdf("Density_CD8_Stress.pdf",width=5,height=3)
p = ggplot(dat, aes(x=dpt_order_perc, fill=meta.cluster)) +
      geom_density(size=0.2, adjust=1.5, n=50, trim=F, color="black", alpha=0.75) +
      theme_classic2() +
      scale_fill_manual(values=colSet[[1]][unique(dat$meta.cluster)])
print(p)
dens_p = p + theme_void() + theme(plot.margin=margin(t=0,r=0,b=0,l=0,unit="cm"),legend.position="none")+xlim(0,1)+scale_y_continuous(expand = c(0,0))+scale_x_continuous(expand = c(0,0))
dev.off()

### DEG
assay(sce, "exprs") <- CD8@assays$RNA@scale.data[rownames(assay(sce,"logcounts")),colnames(assay(sce,"logcounts"))]
gam.test = testProcess(sce, tfs)
saveRDS(gam.test, file="CD8.dyn_genes.rds")

genes <- intersect(tfs, rownames(assay(sce,"exprs")))
a <- assay(sce, "exprs")[genes,]
### heatmap
p.cutoff = 0.05
c.cutoff = 0.2
#
gam.test$sig = ifelse(gam.test$adj.p < p.cutoff & abs(gam.test$coef) > c.cutoff, T, F)
#
plot.gene = c(as.character(gam.test[gam.test$sig, "geneSymbol"]))
highlight.genes = c("GZMB","PRF1","DNAJB1","AHR","HSPA1B","HSPH1","EPAS1","HIF1A")
annot=rowAnnotation(gene=anno_mark(at=match(highlight.genes, plot.gene), labels=highlight.genes))
colSet[['dpt_order_perc']]=colorRamp2(c(0,max(sce$dpt_order_perc)),c("white","darkblue"))
#
set.seed(1)
sce$meta.cluster = as.character(sce$meta.cluster)
sce$cluster.name = as.character(sce$cluster.name)
sce$cluster.name = factor(sce$cluster.name)
pdf("CD8_All_gene_Stress.pdf",width=7,height=8)
hm.obj = heatmap_sm(sce[plot.gene,], assay.name="exprs",
                 colSet=colSet,
                 out.prefix="CD8.dyn_genes.heatmap",
                 columns=c("dpt_order_perc"), 
                 ncell.downsample=500,
                 use_raster=T, raster_quality=2,
                 columns.order="dpt_order_perc",
                 show_column_names=F, show_row_names=F,
                 row_names_side="right",
                 row_gap = unit(1, "mm"), column_gap = unit(0, "mm"),
                 border=F, ann.bar.height=0.5, palette.name="blue2green2red",
                 pdf.width=10.5, pdf.height=7, do.scale=F,
                 z.lo=-0.5, z.hi=0.5, z.step=0.5, 
                 do.clustering.row=T, do.clustering.col=F,
                 dend.row=T,  right_annotation=annot,
                 clustering.distance="cosine",
                 clustering.method="ward.D2",
                 #row_split=5
                 row_km=3, row_km_repeats=200, #row_title=NULL
                 )
print(hm.obj)
dev.off()

pdf("CD8_genes.heatmap_addDensity_All_gene_Stress.pdf", width=5.5, height=4)
HeatmapAnnotation(ggplot=anno_empty(height=unit(1.2, "cm")), which="column") %v% hm.obj

decorate_annotation("ggplot", {
  vp = current.viewport()$name
  print(dens_p, vp=vp)
})
dev.off()


#CD4
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/Tcell/Glioma/trajectory/CD4")

seurat <- readRDS("../../GBM_Tcell_anno.rds")
CD4 <- subset(seurat,idents=c("CD4_Naive","CD4_Memory","CD4_Stress"))
DefaultAssay(CD4) <- "RNA"
CD4 <- ScaleData(CD4, features = rownames(CD4))
sce = as.SingleCellExperiment(CD4)

diffmap <- CD4@meta.data
rownames(diffmap) <- diffmap$barcode
diffmap <- diffmap[,c("UMAP_1","UMAP_2")]

dpt = read.csv("CD4_Slingshot_pseudotime.xls",sep="\t")
diffmap$dpt_pseudotime = dpt[colnames(CD4),1]

sce = sce[,rownames(diffmap)]
reducedDim(sce, "diffmap") = diffmap[,c("UMAP_1","UMAP_2")]
sce$dpt_pseudotime = diffmap$dpt_pseudotime
sce$meta.cluster = as.character(sce$anno)
#
sce = convertDPTperc(sce)
sce$cluster.name = sce$anno

####
colSet = list()
color = c("#BC3C29FF", "#ADB17DFF", "#B1746FFF")
names(color) <- c("CD4_Naive","CD4_Memory","CD4_Stress")
colSet[[1]] <- color

### density
dat = colData(sce)[,c("dpt_pseudotime","dpt_order_perc","meta.cluster")] %>% as.data.frame
pdf("Density.pdf",width=5,height=3)
p = ggplot(dat, aes(x=dpt_order_perc, fill=meta.cluster)) +
      geom_density(size=0.2, adjust=1.5, n=50, trim=F, color="black", alpha=0.75) +
      theme_classic2() +
      scale_fill_manual(values=colSet[[1]][unique(dat$meta.cluster)])
print(p)
dens_p = p + theme_void() + theme(plot.margin=margin(t=0,r=0,b=0,l=0,unit="cm"),legend.position="none")+xlim(0,1)+scale_y_continuous(expand = c(0,0))+scale_x_continuous(expand = c(0,0))
dev.off()

### DEG
assay(sce, "exprs") <- CD4@assays$RNA@scale.data[rownames(assay(sce,"logcounts")),colnames(assay(sce,"logcounts"))]
gam.test = testProcess(sce, tfs)
saveRDS(gam.test, file="CD4.dyn_genes.rds")

genes <- intersect(tfs, rownames(assay(sce,"exprs")))
a <- assay(sce, "exprs")[genes,]

### heatmap
p.cutoff = 0.05
c.cutoff = 0.2
#
gam.test$sig = ifelse(gam.test$adj.p < p.cutoff & abs(gam.test$coef) > c.cutoff, T, F)
#
plot.gene = as.character(gam.test[gam.test$sig, "geneSymbol"])
highlight.genes = c("LEF1","CCR7","TCF7","SELL","KLRB1","RORA","GPR183","S100A4","IL2RB","HSPB1","HSPA9","HSPA8","HSPA6","HSPD1","HSP90AB1")
annot=rowAnnotation(gene=anno_mark(at=match(highlight.genes, plot.gene), labels=highlight.genes))
colSet[['dpt_order_perc']]=colorRamp2(c(0,max(sce$dpt_order_perc)),c("white","darkblue"))
#colSet[['cluster.name']] = structure(nam.conv$mcolor, names=sce$cluster.name)
#
set.seed(1)
sce$meta.cluster = as.character(sce$meta.cluster)
sce$cluster.name = as.character(sce$cluster.name)
sce$cluster.name = factor(sce$cluster.name)
pdf("CD4_All_gene.pdf",width=8,height=8)
hm.obj = heatmap_sm(sce[plot.gene,], assay.name="exprs",
                 colSet=colSet,
                 out.prefix="CD4.dyn_genes.heatmap",
                 columns=c("dpt_order_perc"), 
                 ncell.downsample=500,
                 use_raster=T, raster_quality=2,
                 columns.order="dpt_order_perc",
                 show_column_names=F, show_row_names=F,
                 row_names_side="right",
                 row_gap = unit(1, "mm"), column_gap = unit(0, "mm"),
                 border=F, ann.bar.height=0.5, palette.name="blue2green2red",
                 pdf.width=10.5, pdf.height=7, do.scale=F,
                 z.lo=-0.5, z.hi=0.5, z.step=0.5, 
                 do.clustering.row=T, do.clustering.col=F,
                 dend.row=T,  right_annotation=annot,
                 clustering.distance="cosine",
                 clustering.method="ward.D2",
                 #row_split=5
                 row_km=3, row_km_repeats=200, #row_title=NULL
                 )
print(hm.obj)
dev.off()


pdf("CD4_genes.heatmap_addDensity_All_gene_V2.pdf", width=5.5, height=4)
HeatmapAnnotation(ggplot=anno_empty(height=unit(1.2, "cm")), which="column") %v% hm.obj

decorate_annotation("ggplot", {
  vp = current.viewport()$name
  print(dens_p, vp=vp)
})
dev.off()
