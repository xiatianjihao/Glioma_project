library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/spatial_data_2023/NC/GSE194329/result")
dir <- '/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/spatial_data_2023/NC/GSE194329'
#expr.url <- paste0(dir,"/UKF242_T_ST/outs/filtered_feature_bc_matrix")
#expr.data <- Seurat::Read10X(expr.url)
#UKF242_T_ST <- Seurat::CreateSeuratObject(counts = expr.data, project = 'UKF242_T_ST', assay = 'Spatial')
#UKF242_T_ST$slice <- 1
#UKF242_T_ST$region <- 'UKF242_T_ST'
# Load the image data
#img.url <- paste0(dir,"/UKF242_T_ST/outs/spatial")
#img <- Seurat::Read10X_Image(image.dir = img.url)
#Seurat::DefaultAssay(object = img) <- 'Spatial'
#img <- img[colnames(x = UKF242_T_ST)]
#UKF242_T_ST[['image']] <- img

files <- list.files(path=dir,pattern="spaceranger_out")


for( id in files){
	dat_dir <- paste0(dir,"/",id)
	object <- Seurat::Load10X_Spatial(dat_dir)
	plot1 <- VlnPlot(object, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
	plot2 <- SpatialFeaturePlot(object, features = "nCount_Spatial") + theme(legend.position = "right")
	ggsave(wrap_plots(plot1, plot2),file=paste0(id,"_QC.pdf"),width=10,height=5)
	object <- SCTransform(object, assay = "Spatial",return.only.var.genes = FALSE, verbose = FALSE)
	saveRDS(object,file=paste0(id,".rds"))
	p1 <- SpatialFeaturePlot(object, features = c("CX3CR1","TMEM119","P2RY12","BIN1","NAV3","SPP1","CD3D","CD8A","CD8B","CD44", "ITGA4"),slot="data", image.alpha=0.5)
	ggsave(p1, file=paste0(id,"_SpatialFeaturePlot.data.pdf"),width=16,height=16)
	p2 <- SpatialFeaturePlot(object, features = c("CX3CR1","TMEM119","P2RY12","BIN1","NAV3","SPP1","CD3D","CD8A","CD8B","CD44", "ITGA4"),slot="scale.data", image.alpha=0.5)
	ggsave(p2, file=paste0(id,"_SpatialFeaturePlot.scaledata.pdf"),width=16,height=16)
	expr <- t(as.matrix(object@assays$Spatial@data))

	all_num <- length(rownames(expr))

	mg <- expr[which(expr[,"CX3CR1"] >0 | expr[,"TMEM119"] >0 | expr[,"P2RY12"] >0),]

	mg_ratio<- length(rownames(mg))/all_num

	mg_vs_spp1 <- mg[which(mg[,"SPP1"] > 0),]
	mg_vs_spp1_ratio <- length(rownames(mg_vs_spp1))/all_num

	CD8 <- expr[which(expr[,"CD8A"] >0 | expr[,"CD8B"] >0),]
	CD8_ratio <- length(rownames(CD8))/all_num

	CD8_vs_ITGA4 <- CD8[which(CD8[,"ITGA4"] > 0),]
	CD8_vs_ITGA4_ratio <- length(rownames(CD8[which(CD8[,"ITGA4"] > 0),]))/all_num

	CD8_vs_CD44 <- CD8[which(CD8[,"CD44"] > 0),]
	CD8_vs_CD44_ratio <- length(rownames(CD8[which(CD8[,"CD44"] > 0),]))/all_num


	mg_spp1_vs_CD8_ITGA4_ratio <- length(rownames(mg_vs_spp1[intersect(rownames(mg_vs_spp1),rownames(CD8_vs_ITGA4)),c("CX3CR1","TMEM119","P2RY12","SPP1","CD8A","CD8B","ITGA4","CD44")]))/all_num
	mg_spp1_vs_CD8_CD44_ratio <- length(rownames(mg_vs_spp1[intersect(rownames(mg_vs_spp1),rownames(CD8_vs_CD44)),c("CX3CR1","TMEM119","P2RY12","SPP1","CD8A","CD8B","ITGA4","CD44")]))/all_num



	df <- data.frame(Group=c("MG","MG-SPP1+","CD8","CD8-ITGA4+","CD8-CD44+","MG-SPP1+CD8-ITGA4+","MG-SPP1+CD8-CD44+"),Fraction=c(mg_ratio,mg_vs_spp1_ratio,CD8_ratio,CD8_vs_ITGA4_ratio,CD8_vs_CD44_ratio,mg_spp1_vs_CD8_ITGA4_ratio,mg_spp1_vs_CD8_CD44_ratio))
	p <- ggplot(df,aes(x=Group,y=Fraction))+geom_bar(stat="identity",fill="lightblue")+geom_text(aes(label=sprintf("%.3f",Fraction)),vjust=-0.2)+theme_classic()+theme(axis.text = element_text(size=12,color="black"),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 12),legend.text = element_text(size = 12))
	ggsave(p, file=paste0(id,"_spots.pdf"),width=6,height=6)
}


#merge all samples
GBM1 <- readRDS("GBM1_spaceranger_out.rds")
GBM2 <- readRDS("GBM2_spaceranger_out.rds")
GBM3 <- readRDS("GBM3_spaceranger_out.rds")
GBM4 <- readRDS("GBM4_spaceranger_out.rds")
GBM5_1 <- readRDS("GBM5_1_spaceranger_out.rds")
GBM5_2 <- readRDS("GBM5_2_spaceranger_out.rds")

#60124 features across 16344 samples within 2 assays
brain.merge <- merge(GBM1,y=c(GBM2,GBM3,GBM4,GBM5_1,GBM5_2), add.cell.ids=c("GBM1","GBM2","GBM3","GBM4","GBM5_1","GBM5_2"))

DefaultAssay(brain.merge) <- "SCT"

saveRDS(brain.merge,file="brain.merge.rds")

p1 <- SpatialFeaturePlot(brain.merge, features = c("CX3CR1","TMEM119","P2RY12","SPP1","CD8A","CD8B","CD44","ITGA4"),slot="data", image.alpha=0.5)
ggsave(p1, file="merge_SpatialFeaturePlot.data.pdf",width=20,height=20)

expr <- t(as.matrix(brain.merge@assays$Spatial@data))

all_num <- length(rownames(expr))

mg <- expr[which(expr[,"CX3CR1"] >0 | expr[,"TMEM119"] >0 | expr[,"P2RY12"] >0),]

mg_ratio<- length(rownames(mg))/all_num

mg_vs_spp1 <- mg[which(mg[,"SPP1"] > 0),]
mg_vs_spp1_ratio <- length(rownames(mg_vs_spp1))/all_num

CD8 <- expr[which(expr[,"CD8A"] >0 | expr[,"CD8B"] >0),]
CD8_ratio <- length(rownames(CD8))/all_num

CD8_vs_ITGA4 <- CD8[which(CD8[,"ITGA4"] > 0),]
CD8_vs_ITGA4_ratio <- length(rownames(CD8[which(CD8[,"ITGA4"] > 0),]))/all_num

CD8_vs_CD44 <- CD8[which(CD8[,"CD44"] > 0),]
CD8_vs_CD44_ratio <- length(rownames(CD8[which(CD8[,"CD44"] > 0),]))/all_num


mg_spp1_vs_CD8_ITGA4_ratio <- length(rownames(mg_vs_spp1[intersect(rownames(mg_vs_spp1),rownames(CD8_vs_ITGA4)),c("CX3CR1","TMEM119","P2RY12","SPP1","CD8A","CD8B","ITGA4","CD44")]))/all_num
mg_spp1_vs_CD8_CD44_ratio <- length(rownames(mg_vs_spp1[intersect(rownames(mg_vs_spp1),rownames(CD8_vs_CD44)),c("CX3CR1","TMEM119","P2RY12","SPP1","CD8A","CD8B","ITGA4","CD44")]))/all_num



df <- data.frame(Group=c("MG","MG-SPP1+","CD8","CD8-ITGA4+","CD8-CD44+","MG-SPP1+CD8-ITGA4+","MG-SPP1+CD8-CD44+"),Fraction=c(mg_ratio,mg_vs_spp1_ratio,CD8_ratio,CD8_vs_ITGA4_ratio,CD8_vs_CD44_ratio,mg_spp1_vs_CD8_ITGA4_ratio,mg_spp1_vs_CD8_CD44_ratio))
p <- ggplot(df,aes(x=Group,y=Fraction))+geom_bar(stat="identity",fill="lightblue")+geom_text(aes(label=sprintf("%.3f",Fraction)),vjust=-0.2)+theme_classic()+theme(axis.text = element_text(size=12,color="black"),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 12),legend.text = element_text(size = 12))
ggsave(p, file="MG_merge_spots.pdf",width=6,height=6)


#########MDM SPP1 vs CD44/ITGA4
brain.merge <- readRDS("brain.merge.rds")


expr <- t(as.matrix(brain.merge@assays$Spatial@data))      #16344 samplesß

all_num <- length(rownames(expr))

MDM <- expr[which(expr[,"ITGA4"] >0 | expr[,"MRC1"] >0 | expr[,"TGFB1"] >0 | expr[,"THBD"] > 0),]

MDM_ratio<- length(rownames(MDM))/all_num

MDM_vs_spp1 <- MDM[which(MDM[,"SPP1"] > 0),]
MDM_vs_spp1_ratio <- length(rownames(MDM_vs_spp1))/all_num

CD8 <- expr[which(expr[,"CD8A"] >0 | expr[,"CD8B"] >0),]
CD8_ratio <- length(rownames(CD8))/all_num

CD8_vs_ITGA4 <- CD8[which(CD8[,"ITGA4"] > 0),]
CD8_vs_ITGA4_ratio <- length(rownames(CD8[which(CD8[,"ITGA4"] > 0),]))/all_num

CD8_vs_CD44 <- CD8[which(CD8[,"CD44"] > 0),]
CD8_vs_CD44_ratio <- length(rownames(CD8[which(CD8[,"CD44"] > 0),]))/all_num


MDM_spp1_vs_CD8_ITGA4_ratio <- length(rownames(MDM_vs_spp1[intersect(rownames(MDM_vs_spp1),rownames(CD8_vs_ITGA4)),c("ITGA4","MRC1","TGFB1","THBD","SPP1","CD8A","CD8B","ITGA4","CD44")]))/all_num
MDM_spp1_vs_CD8_CD44_ratio <- length(rownames(MDM_vs_spp1[intersect(rownames(MDM_vs_spp1),rownames(CD8_vs_CD44)),c("ITGA4","MRC1","TGFB1","THBD","SPP1","CD8A","CD8B","ITGA4","CD44")]))/all_num



df <- data.frame(Group=c("MDM","MDM-SPP1+","CD8","CD8-ITGA4+","CD8-CD44+","MDM-SPP1+CD8-ITGA4+","MDM-SPP1+CD8-CD44+"),Fraction=c(MDM_ratio,MDM_vs_spp1_ratio,CD8_ratio,CD8_vs_ITGA4_ratio,CD8_vs_CD44_ratio,MDM_spp1_vs_CD8_ITGA4_ratio,MDM_spp1_vs_CD8_CD44_ratio))
p <- ggplot(df,aes(x=Group,y=Fraction))+geom_bar(stat="identity",fill="lightblue")+geom_text(aes(label=sprintf("%.3f",Fraction)),vjust=-0.2)+theme_classic()+theme(axis.text = element_text(size=12,color="black"),axis.text.x = element_text(size = 12,angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 12),legend.text = element_text(size = 12))
ggsave(p, file="MDM_merge_spots.pdf",width=6,height=6)

