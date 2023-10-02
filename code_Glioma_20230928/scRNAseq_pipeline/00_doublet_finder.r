library(Seurat)
library(DoubletFinder)
library(ggplot2)

dirs <- list.files(path="/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/data",pattern="*")
samples <- c("GBM1","PBMC1","GBM2","PBMC2","GBM3","PBMC3","GBM4","PBMC4","GBM5","PBMC5","GBM6","PBMC6","GBM7","PBMC7","GBM8","PBMC8","GBM9","PBMC9")

source("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/code_Glioma_project_upload/scRNAseq_pipeline/comm_func.r")


sceList <- lapply(dirs, function(x){
	name <- samples[which(dirs==x)]
	print(x)
	print(name)
	dir.create(paste0("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/result/",name))
	input_dir <- paste0("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/data/",x,"/outs/filtered_feature_bc_matrix")
	out_dir <- paste0("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/result/",name)
	QC_for_singlesample(input_dir,out_dir,name)
})

sceList <- lapply(samples, function(x){
	out_dir <- paste0("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/result/",x)
	DoubletFinder_for_singlesample(out_dir,x)
})

