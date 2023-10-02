
SoupX_DoubletFinder_for_singlesample <- function(raw_dir,filter_dir,out_dir,name){
	raw <- Read10X_h5(raw_dir)
	filter <- Read10X(filter_dir)
	srat <- CreateSeuratObject(counts =filter)
	soup.channel <- SoupChannel(raw,filter)
	srat <- SCTransform(srat,verbose=F)
	srat <- RunPCA(srat)
	srat <- RunUMAP(srat,dims=1:30)
	srat <- FindNeighbors(srat,dims=1:30)
	srat <- FindClusters(srat)
	meta <- srat@meta.data
	umap <- srat@reductions$umap@cell.embeddings
	soup.channel <-setClusters(soup.channel,setNames(meta$seurat_clusters,rownames(meta)))
	soup.channel <- setDR(soup.channel,umap)
	soup.channel <- autoEstCont(soup.channel)
	adj.matrix <- adjustCounts(soup.channel,roundToInt=T)

	DropletUtils::write10xCounts(name,adj.matrix)
	data <- CreateSeuratObject(counts = adj.matrix, project = name, min.cells = 3, min.features = 200)
	data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
	ribo_genes <- rownames(data@assays$RNA@data)[grep("^RP[SL]",rownames(data@assays$RNA@data))]
	data[["percent.ribo"]] <- PercentageFeatureSet(data, features=ribo_genes)
	hb_genes <- rownames(data@assays$RNA@data)[grep("^HB[^(P)]",rownames(data@assays$RNA@data))]
	data[["percent.hb"]] <- PercentageFeatureSet(data, features=hb_genes)
	p1 <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo","percent.hb"), ncol = 5, pt.size=0.01)
	ggsave(p1, file=paste0(out_dir,"/",name,"_QC.pdf"),width=15,height=4)
	p1 <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo","percent.hb"), ncol = 5, pt.size=0)
	ggsave(p1, file=paste0(out_dir,"/",name,"_QC_v2.pdf"),width=15,height=4)
	saveRDS(data,file=paste0(out_dir,"/",name,"_all.rds"))

	#data <- readRDS(paste0(out_dir,"/",name,"_all.rds"))
	data <- subset(data, subset=nFeature_RNA>200 & nFeature_RNA < 6500 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mt < 15)
	data <- NormalizeData(data)
	data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
	data <- ScaleData(data)
	data <- RunPCA(data)
	data <- RunUMAP(data, dims = 1:20)
	sweep.res.list_data <- paramSweep_v3(data, PCs = 1:20, sct = FALSE)
	sweep.stats_data <- summarizeSweep(sweep.res.list_data, GT = FALSE)
	pdf(paste0(out_dir,"/",name,"_FindPK.pdf"))
	bcmvn_data <- find.pK(sweep.stats_data)
	dev.off()
	mpK<-as.numeric(as.vector(bcmvn_data$pK[which.max(bcmvn_data$BCmetric)]))

	cellnum <- length(rownames(data@meta.data))
	if(cellnum < 5000){
		doublet_rate <- 0.075
	}else if(cellnum > 5000 && cellnum < 10000){
		doublet_rate  <- 0.10
	}else if(cellnum >10000 && cellnum < 20000){
		doublet_rate <- 0.15
	}else{
		doublet_rate <- 0.20
	}
	## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
	nExp_poi <- round(doublet_rate *length(rownames(data@meta.data)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset

	data <- doubletFinder_v3(data, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	colnames(data@meta.data) <- c("orig.ident","nCount_RNA", "nFeature_RNA","percent.mt", "percent.ribo","percent.hb","pANN","DN")
	saveRDS(data, file=paste0(out_dir,"/",name,"_finddoublets.rds"))

	data <- SetIdent(data,value="DN")
	p1 <- DimPlot(data, reduction = "umap", cols = c("red","black"))+theme(text = element_text(size=12),axis.text.x = element_text(color="black"))+xlab("UMAP-1")+ylab("UMAP-2")
	ggsave(p1, file=paste0(out_dir,"/",name,"_data_DoubletFinder.pdf"),width=5,height=4)

	p2 <- FeaturePlot(object = data, features=c("pANN"))
	ggsave(p2, file=paste0(out_dir,"/",name,'_DN_value.pdf'),width=5,height=4)


	data2 <- subset(data,ident="Singlet")
	saveRDS(data2,file=paste0(out_dir,"/",name,"_singlet.rds"))
}


## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = unit(c(-0.10, 0, -0.10, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 0),
          axis.text.y = element_text(size = rel(1)),
          plot.margin = plot.margin )
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.10, 0, -0.10, 0), "cm"),
                          ...) {

  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))

  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.ticks.x = element_line())

  # change the y-axis tick to only max value
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
                            scale_y_continuous(breaks = c(0, y/2, y)) +
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

