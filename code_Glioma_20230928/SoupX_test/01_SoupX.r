library(Seurat)
library(SoupX)
library(DropletUtils)
library(DoubletFinder)
library(knitr)
library(ggplot2)


dirs <- list.files(path="/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/data",pattern="[0-9]*-PBMC|GBM")
samples <- c("GBM1","PBMC1","GBM2","PBMC2","GBM3","PBMC3","GBM4","PBMC4","GBM5","PBMC5","GBM6","PBMC6","GBM7","PBMC7","GBM8","PBMC8","GBM9","PBMC9")

source("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/20230919_revision_code/comm_func.r")


sceList <- lapply(dirs, function(x){
	name <- samples[which(dirs==x)]
	print(x)
	print(name)
	raw_dir <- paste0("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/data/",x,"/outs/raw_feature_bc_matrix.h5")
	filter_dir <- paste0("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/data/",x,"/outs/filtered_feature_bc_matrix")
	out_dir <- paste0("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/20230919_revision_code/01_SoupX/",name,"_out")
	dir.create(out_dir)
	SoupX_DoubletFinder_for_singlesample(raw_dir,filter_dir,out_dir,name)
})

