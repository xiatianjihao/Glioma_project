library(tidyverse)
library(RColorBrewer)
library(scales)
library(igraph)
library(stringr)


#PBMC network
pvalues=read.table("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/cellphonedb/V3/compare_GBM_vs_PBMC/PBMC_pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.05))
colnames(statdf)=c("number")

statdf$indexb=str_replace(rownames(statdf),"^.*\\.","")
statdf$indexa=str_replace(rownames(statdf),"\\..*$","")
rankname=sort(unique(statdf$indexa))

A=c()
B=c()
C=c()
remaining=rankname
for (i in rankname[-9]) {
  remaining=setdiff(remaining,i)
  for (j in remaining) {
    count=statdf[statdf$indexa == i & statdf$indexb == j,"number"]+
      statdf[statdf$indexb == i & statdf$indexa == j,"number"]
    A=append(A,i)
    B=append(B,j)
    C=append(C,count)
  }
}

statdf2=data.frame(indexa=A,indexb=B,number=C)
statdf2=statdf2 %>% rbind(statdf[statdf$indexa==statdf$indexb,c("indexa","indexb","number")])
statdf2=statdf2[statdf2$number > 0,] #过滤掉值为0的观测

#设置节点和连线的颜色
color1=c("#4DBBD5B2","#FDB462", "#B3DE69","#B09C85B2", "#BC80BD","#F39B7FB2","#00A087B2")
names(color1)=rankname
color2=colorRampPalette(brewer.pal(9, "Reds")[3:7])(20) #将颜色分成多少份，取决于互作关系数目的最大值
names(color2)=1:20 #每一份颜色用对应的数字命名

#做网络图
##下面的四行代码相对固定
net <- graph_from_data_frame(statdf2[,c("indexa","indexb","number")])
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = order(membership(group)))

E(net)$width <- E(net)$number / 50 #将数值映射到连线的宽度，有时还需要微调，这里除以2就是这个目的
E(net)$color <- color2[as.character(ifelse(E(net)$number > 20,20,E(net)$number))] #用前面设置好的颜色赋给连线，颜色深浅对应数值大小
E(net)$label = E(net)$number #连线的标注
E(net)$label.color <- "black" #连线标注的颜色
V(net)$label.color <- "black" #节点标注的颜色
V(net)$color <- color1[names(V(net))] #节点的填充颜色，前面已经设置了；V(net)返回节点信息

#调整节点位置的线条角度
##如果没有这两行代码，节点位置的圆圈是向右的
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]

pdf("interaction.num.PBMC.pdf",width = 6,height = 6)
plot(net,
     edge.arrow.size = 1, #连线不带箭头
     edge.curved = 0, #连线不弯曲
     vertex.frame.color = "black", #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = 30) #节点大小
dev.off()

#Glioma network
pvalues=read.table("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/cellphonedb/V3/compare_GBM_vs_PBMC/GBM_pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.05))
colnames(statdf)=c("number")

statdf$indexb=str_replace(rownames(statdf),"^.*\\.","")
statdf$indexa=str_replace(rownames(statdf),"\\..*$","")
rankname=sort(unique(statdf$indexa))

A=c()
B=c()
C=c()
remaining=rankname
for (i in rankname[-8]) {
  remaining=setdiff(remaining,i)
  for (j in remaining) {
    count=statdf[statdf$indexa == i & statdf$indexb == j,"number"]+
      statdf[statdf$indexb == i & statdf$indexa == j,"number"]
    A=append(A,i)
    B=append(B,j)
    C=append(C,count)
  }
}

statdf2=data.frame(indexa=A,indexb=B,number=C)
statdf2=statdf2 %>% rbind(statdf[statdf$indexa==statdf$indexb,c("indexa","indexb","number")])
statdf2=statdf2[statdf2$number > 0,] #过滤掉值为0的观测

#"B"        "CD4"      "CD8"      "DC"       "MDM"      "MG"       "Monocyte" "NK"
#"#4DBBD5B2","#FDB462", "#B3DE69","#B09C85B2", "#BC80BD","#D9D9D9","#F39B7FB2","#00A087B2"

#设置节点和连线的颜色
#color1=c("#8DD3C7", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD","#4B0082","#00FFFF","#FFA500")
color1=c("#4DBBD5B2","#FDB462", "#B3DE69","#B09C85B2", "#BC80BD","#D9D9D9","#F39B7FB2","#00A087B2")
names(color1)=rankname
color2=colorRampPalette(brewer.pal(9, "Reds")[3:7])(20) #将颜色分成多少份，取决于互作关系数目的最大值
names(color2)=1:20 #每一份颜色用对应的数字命名

#做网络图
##下面的四行代码相对固定
net <- graph_from_data_frame(statdf2[,c("indexa","indexb","number")])
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = order(membership(group)))

E(net)$width <- E(net)$number / 50 #将数值映射到连线的宽度，有时还需要微调，这里除以2就是这个目的
E(net)$color <- color2[as.character(ifelse(E(net)$number > 20,20,E(net)$number))] #用前面设置好的颜色赋给连线，颜色深浅对应数值大小
E(net)$label = E(net)$number #连线的标注
E(net)$label.color <- "black" #连线标注的颜色
V(net)$label.color <- "black" #节点标注的颜色
V(net)$color <- color1[names(V(net))] #节点的填充颜色，前面已经设置了；V(net)返回节点信息

#调整节点位置的线条角度
##如果没有这两行代码，节点位置的圆圈是向右的
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]

pdf("interaction.num.Glioma.pdf", width = 5,height = 5)
plot(net,
     edge.arrow.size = 1, #连线不带箭头
     edge.curved = 0, #连线不弯曲
     vertex.frame.color = "black", #节点外框颜色
     layout = coords,
     vertex.label.cex = 1, #节点标注字体大小
     vertex.size = 30) #节点大小
dev.off()


#heatmap
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/cellphonedb/V1/compare_GBM_vs_PBMC")

library(psych)
library(tidyverse)
library(ggtext)
library(glue)
library(RColorBrewer)
library(scales)
library(stringr)

#statics figure
GBM_pvalues=read.table("GBM_pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
GBM_pvalues=GBM_pvalues[,12:dim(GBM_pvalues)[2]] #此时不关注前11列
GBMstatdf=as.data.frame(colSums(GBM_pvalues < 0.05)) #统计在某一种细胞pair的情况之下，显著的受配体pair的数目；阈值可以自己选
colnames(GBMstatdf)=c("number")
rownames(GBMstatdf)=gsub("MDM\\.","MDM-",rownames(GBMstatdf))
rownames(GBMstatdf)=gsub("MG\\.","MG-",rownames(GBMstatdf))

#排在前面的分子定义为indexa；排在后面的分子定义为indexb
GBMstatdf$indexb=str_replace(rownames(GBMstatdf),"^.*\\.","")
GBMstatdf$indexa=str_replace(rownames(GBMstatdf),"\\..*$","")
#设置合适的细胞类型的顺序
rankname=sort(unique(GBMstatdf$indexa)) 
#转成因子类型，画图时，图形将按照预先设置的顺序排列
GBMstatdf$indexa=factor(GBMstatdf$indexa,levels = rankname)
GBMstatdf$indexb=factor(GBMstatdf$indexb,levels = rankname)

GBMstatdf %>% ggplot(aes(x=indexa,y=indexb,fill=number))+geom_tile(color="white")+
  geom_text(aes(label=number))+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits = c(0,320))+
  scale_x_discrete("cluster 1 produces molecule 1")+
  scale_y_discrete("cluster 2 produces molecule 2")+
  theme_minimal()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = NULL, angle = 45, size=12, color="black"),
    axis.text.y = element_text(size=12, color="black"),
    panel.grid = element_blank()
  )
ggsave(filename = "Glioma.interaction.num.pdf",device = "pdf",width = 8,height = 6)

GBMstatdf2 <- GBMstatdf[grep("^CD",rownames(GBMstatdf)),]
GBMstatdf3 <- GBMstatdf2[-grep("B$|NK$|Treg$|Stress$|Naive$|Memory$|Cytotoxicity$|Exhaustion$",rownames(GBMstatdf2)),]

GBMstatdf3 %>% ggplot(aes(x=indexa,y=indexb,fill=number))+geom_tile(color="white")+
  geom_text(aes(label=number))+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits = c(50,150))+
  scale_x_discrete("cluster 1 produces molecule 1")+
  scale_y_discrete("cluster 2 produces molecule 2")+
  theme_minimal()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = NULL, angle = 45, size=12, color="black"),
    axis.text.y = element_text(size=12, color="black"),
    panel.grid = element_blank()
  )
ggsave(filename = "Glioma.interaction.num_v2.pdf",device = "pdf",width = 4,height = 3)



#statics figure
PBMC_pvalues=read.table("PBMC_pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
PBMC_pvalues=PBMC_pvalues[,12:dim(PBMC_pvalues)[2]] #此时不关注前11列
PBMCstatdf=as.data.frame(colSums(PBMC_pvalues < 0.05)) #统计在某一种细胞pair的情况之下，显著的受配体pair的数目；阈值可以自己选
colnames(PBMCstatdf)=c("number")

#排在前面的分子定义为indexa；排在后面的分子定义为indexb
PBMCstatdf$indexb=str_replace(rownames(PBMCstatdf),"^.*\\.","")
PBMCstatdf$indexa=str_replace(rownames(PBMCstatdf),"\\..*$","")
#设置合适的细胞类型的顺序
rankname=sort(unique(PBMCstatdf$indexa)) 
#转成因子类型，画图时，图形将按照预先设置的顺序排列
PBMCstatdf$indexa=factor(PBMCstatdf$indexa,levels = rankname)
PBMCstatdf$indexb=factor(PBMCstatdf$indexb,levels = rankname)


PBMCstatdf %>% ggplot(aes(x=indexa,y=indexb,fill=number))+geom_tile(color="white")+
  geom_text(aes(label=number))+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits = c(0,320))+
  scale_x_discrete("cluster 1 produces molecule 1")+
  scale_y_discrete("cluster 2 produces molecule 2")+
  theme_minimal()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = NULL, angle = 45, size=12, color="black"),
    axis.text.y = element_text(size=12, color="black"),
    panel.grid = element_blank()
  )
ggsave(filename = "PBMC.interaction.num.pdf",device = "pdf",width = 8,height = 6)


PBMCstatdf2 <- PBMCstatdf[grep("^CD",rownames(PBMCstatdf)),]
PBMCstatdf3 <- PBMCstatdf2[-grep("B$|NK$|Treg$|Stress$|Naive$|memory$|Cytotoxicity$",rownames(PBMCstatdf2)),]

PBMCstatdf3 %>% ggplot(aes(x=indexa,y=indexb,fill=number))+geom_tile(color="white")+
  geom_text(aes(label=number))+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits = c(50,150))+
  scale_x_discrete("cluster 1 produces molecule 1")+
  scale_y_discrete("cluster 2 produces molecule 2")+
  theme_minimal()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = NULL, angle = 45, size=12, color="black"),
    axis.text.y = element_text(size=12, color="black"),
    panel.grid = element_blank()
  )
ggsave(filename = "PBMC.interaction.num_v2.pdf",device = "pdf",width = 4,height = 2.35)


#CD8_Cytotoxicity vs myeloid interaction
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/cellphonedb/V1/compare_GBM_vs_PBMC")

library(psych)
library(tidyverse)
library(ggtext)
library(glue)
library(RColorBrewer)
library(scales)

#selection
GBMmeans <- read.csv("GBM_means.txt",header=T,sep="\t")
GBMpvals <- read.csv("GBM_pvalues.txt",header=T,sep="\t")

GBMmeans %>% dplyr::select("interacting_pair",starts_with(c("Monocyte","MDM.1","MDM.2","MDM.3","MG.1","MG.2","MG.3","DC")))  %>% reshape2::melt() -> GBMmeansdf

colnames(GBMmeansdf)<- c("interacting_pair","CC","means")

GBMpvals %>% dplyr::select("interacting_pair",starts_with(c("Monocyte","MDM.1","MDM.2","MDM.3","MG.1","MG.2","MG.3","DC")))%>% reshape2::melt()-> GBMpvalsdf

colnames(GBMpvalsdf)<- c("interacting_pair","CC","pvals")

GBMpvalsdf$joinlab<- paste0(GBMpvalsdf$interacting_pair,"_",GBMpvalsdf$CC)
GBMmeansdf$joinlab<- paste0(GBMmeansdf$interacting_pair,"_",GBMmeansdf$CC)
GBMpldf <- merge(GBMpvalsdf,GBMmeansdf,by = "joinlab")

#GBMpldf %>% filter(means > 0) %>%  filter(pvals < 0.05 ) -> filter_result

GBMpldf2 <- GBMpldf[grep("CD8_Cytotoxicity$",GBMpldf$CC.x),]
GBMpldf2$CC.x <- gsub("MDM\\.","MDM-",GBMpldf2$CC.x)
GBMpldf2$CC.x <- gsub("MG\\.","MG-",GBMpldf2$CC.x)

GBMpldf2 %>% filter(pvals < 1) %>% filter(means >1) -> df_v1

#df_v1$CC.x <- paste0("GBM_",df_v1$CC.x)
df_v1$group <- "Glioma"

###################
PBMCmeans <- read.csv("PBMC_means.txt",header=T,sep="\t")
PBMCpvals <- read.csv("PBMC_pvalues.txt",header=T,sep="\t")

PBMCmeans %>% dplyr::select("interacting_pair",starts_with(c("Monocyte","MDM","DC")))  %>% reshape2::melt() -> PBMCmeansdf

colnames(PBMCmeansdf)<- c("interacting_pair","CC","means")

PBMCpvals %>% dplyr::select("interacting_pair",starts_with(c("Monocyte","MDM","DC"))) %>% reshape2::melt()-> PBMCpvalsdf

colnames(PBMCpvalsdf)<- c("interacting_pair","CC","pvals")

PBMCpvalsdf$joinlab<- paste0(PBMCpvalsdf$interacting_pair,"_",PBMCpvalsdf$CC)
PBMCmeansdf$joinlab<- paste0(PBMCmeansdf$interacting_pair,"_",PBMCmeansdf$CC)
PBMCpldf <- merge(PBMCpvalsdf,PBMCmeansdf,by = "joinlab")

#PBMCpldf %>% filter(means > 0) %>%  filter(pvals < 0.05 ) -> filter_result

PBMCpldf2 <- PBMCpldf[grep("CD8_Cytotoxicity$",PBMCpldf$CC.x),]

PBMCpldf2 %>% filter(pvals < 1) %>% filter(means >1) -> df_v2

#df_v1$CC.x <- paste0("GBM_",df_v1$CC.x)
df_v2$group <- "PBMC"

df1 <- rbind(df_v1,df_v2)

df <- df1[which(df1$interacting_pair.x %in% c("HLA-E_KLRD1","HMGB1_CXCR4","TIMP1_CD63","MIF_CXCR4","HLA-E_KLRC1","SPP1_ITGA4","SPP1_CD44","CCL4L2_VSIR","CXCL10_CXCR3")),]

df$CC.x <- paste0(df$group,"_",df$CC.x)

xlab <- c()
i <- 0
for(x in as.vector(df$CC.x)){
  i <- i + 1
  y <- paste0("<span style='color:blue'>",gsub("Glioma_|PBMC_","",str_split(x,"\\.")[[1]][1]),"</span><span style='color:black'> vs </span></span><span style='color:red'>",str_split(x,"\\.")[[1]][2],"</span>")
  xlab[i] <- y
}

ylab <- c()
j <- 0
for(x in as.vector(df$interacting_pair.x)){
  j <- j + 1
  y <- paste0("<span style='color:blue'>",str_split(x,"_")[[1]][1],"</span><span style='color:black'>-</span></span><span style='color:red'>",str_split(x,"_")[[1]][2],"</span>")
  ylab[j] <- y
}

df$xlab <- xlab
df$ylab <- ylab
df$log2means <- log2(df$means)

#df$ylab <- factor(df$ylab,levels=c(df$ylab[c(10,18,26,4,34,38,1,3,42)]))
df$ylab <- factor(df$ylab,levels=c(df$ylab[c(42,3,1,38,34,4,26,18,10)]))

breaks <- c(-2,-1,0,1,2)
pdf("Myeloid_vs_CD8_Cytotoxicity_LR.pdf",width=7.5,height=4)
ggplot(df, aes(x = xlab, y = ylab))+geom_point(aes(size = -log10(pvals+0.0001),color = log2means)) +
geom_vline(aes_string(xintercept = "xlab"),lwd = 0.2,colour = "grey",linetype="dashed") +
geom_hline(aes_string(yintercept = "ylab"),lwd = 0.2,colour = "grey",linetype="dashed")+scale_colour_gradientn(colors = c("#000000","#0000FF","#FFFF00","#FF0000"),breaks = breaks, labels = format(breaks),limits=c(-2,2)) +facet_grid(facets =  ~ group, space = "free", scales = "free")+theme_classic()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1))+theme(axis.text.y = element_markdown(), axis.text.x = element_markdown(angle = 90,hjust = 1,vjust = 0.5)) +
guides(size = guide_legend(title = "-log10(pvalue+1e-04)"),
       color = guide_colorbar(title = "Log2 mean (Molecule1, Molecule2)"))
dev.off()


#Treg_vs_others interaction
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/cellphonedb/V1/compare_GBM_vs_PBMC")

GBMmeans <- read.csv("GBM_means.txt",header=T,sep="\t")
GBMpvals <- read.csv("GBM_pvalues.txt",header=T,sep="\t")

GBMmeans %>% dplyr::select("interacting_pair",starts_with(c("Treg")))  %>% reshape2::melt() -> GBMmeansdf

colnames(GBMmeansdf)<- c("interacting_pair","CC","means")

GBMpvals %>% dplyr::select("interacting_pair",starts_with(c("Treg")))%>% reshape2::melt()-> GBMpvalsdf

colnames(GBMpvalsdf)<- c("interacting_pair","CC","pvals")

GBMpvalsdf$joinlab<- paste0(GBMpvalsdf$interacting_pair,"_",GBMpvalsdf$CC)
GBMmeansdf$joinlab<- paste0(GBMmeansdf$interacting_pair,"_",GBMmeansdf$CC)
GBMpldf <- merge(GBMpvalsdf,GBMmeansdf,by = "joinlab")

GBMpldf$CC.x <- gsub("MDM\\.","MDM-",GBMpldf$CC.x)
GBMpldf$CC.x <- gsub("MG\\.","MG-",GBMpldf$CC.x)

GBMpldf2 <- GBMpldf[-grep("B$|NK$|Treg$",GBMpldf$CC.x),]

GBMpldf2 %>% filter(pvals < 1) %>% filter(means >1) -> df_v1

#df_v1$CC.x <- paste0("GBM_",df_v1$CC.x)
df_v1$group <- "Glioma"

###################
PBMCmeans <- read.csv("PBMC_means.txt",header=T,sep="\t")
PBMCpvals <- read.csv("PBMC_pvalues.txt",header=T,sep="\t")

PBMCmeans %>% dplyr::select("interacting_pair",starts_with(c("Treg")))  %>% reshape2::melt() -> PBMCmeansdf

colnames(PBMCmeansdf)<- c("interacting_pair","CC","means")

PBMCpvals %>% dplyr::select("interacting_pair",starts_with(c("Treg"))) %>% reshape2::melt()-> PBMCpvalsdf

colnames(PBMCpvalsdf)<- c("interacting_pair","CC","pvals")

PBMCpvalsdf$joinlab<- paste0(PBMCpvalsdf$interacting_pair,"_",PBMCpvalsdf$CC)
PBMCmeansdf$joinlab<- paste0(PBMCmeansdf$interacting_pair,"_",PBMCmeansdf$CC)
PBMCpldf <- merge(PBMCpvalsdf,PBMCmeansdf,by = "joinlab")

PBMCpldf2 <- PBMCpldf[-grep("B$|NK$|Treg$",PBMCpldf$CC.x),]

PBMCpldf2 %>% filter(pvals < 1) %>% filter(means >1) -> df_v2

#df_v1$CC.x <- paste0("GBM_",df_v1$CC.x)
df_v2$group <- "PBMC"

df1 <- rbind(df_v1,df_v2)


df1_use <- df1[which(df1$interacting_pair.x %in% c("VIM_CD44","B2M_LILRB2","HLA-E_KLRC1","MIF_CD44")),]


GBMmeans <- read.csv("GBM_means.txt",header=T,sep="\t")
GBMpvals <- read.csv("GBM_pvalues.txt",header=T,sep="\t")

GBMmeans %>% dplyr::select("interacting_pair",starts_with(c("Monocyte","CD4","CD8","MDM","MG","DC")))  %>% reshape2::melt() -> GBMmeansdf

colnames(GBMmeansdf)<- c("interacting_pair","CC","means")

GBMpvals %>% dplyr::select("interacting_pair",starts_with(c("Monocyte","CD4","CD8","MDM","MG","DC")))%>% reshape2::melt()-> GBMpvalsdf

colnames(GBMpvalsdf)<- c("interacting_pair","CC","pvals")

GBMpvalsdf$joinlab<- paste0(GBMpvalsdf$interacting_pair,"_",GBMpvalsdf$CC)
GBMmeansdf$joinlab<- paste0(GBMmeansdf$interacting_pair,"_",GBMmeansdf$CC)
GBMpldf <- merge(GBMpvalsdf,GBMmeansdf,by = "joinlab")

GBMpldf$CC.x <- gsub("MDM\\.","MDM-",GBMpldf$CC.x)
GBMpldf$CC.x <- gsub("MG\\.","MG-",GBMpldf$CC.x)

GBMpldf2 <- GBMpldf[grep("Treg$",GBMpldf$CC.x),]

GBMpldf2 %>% filter(pvals < 1) %>% filter(means >1) -> df_v1

#df_v1$CC.x <- paste0("GBM_",df_v1$CC.x)
df_v1$group <- "Glioma"

###################
PBMCmeans <- read.csv("PBMC_means.txt",header=T,sep="\t")
PBMCpvals <- read.csv("PBMC_pvalues.txt",header=T,sep="\t")

PBMCmeans %>% dplyr::select("interacting_pair",starts_with(c("Monocyte","CD4","CD8","MDM","DC")))  %>% reshape2::melt() -> PBMCmeansdf

colnames(PBMCmeansdf)<- c("interacting_pair","CC","means")

PBMCpvals %>% dplyr::select("interacting_pair",starts_with(c("Monocyte","CD4","CD8","MDM","DC"))) %>% reshape2::melt()-> PBMCpvalsdf

colnames(PBMCpvalsdf)<- c("interacting_pair","CC","pvals")

PBMCpvalsdf$joinlab<- paste0(PBMCpvalsdf$interacting_pair,"_",PBMCpvalsdf$CC)
PBMCmeansdf$joinlab<- paste0(PBMCmeansdf$interacting_pair,"_",PBMCmeansdf$CC)
PBMCpldf <- merge(PBMCpvalsdf,PBMCmeansdf,by = "joinlab")

PBMCpldf2 <- PBMCpldf[grep("Treg$",PBMCpldf$CC.x),]

PBMCpldf2 %>% filter(pvals < 1) %>% filter(means >1) -> df_v2

#df_v1$CC.x <- paste0("GBM_",df_v1$CC.x)
df_v2$group <- "PBMC"

df1 <- rbind(df_v1,df_v2)

df2 <- df1[which(df1$interacting_pair.x %in% c("KLRC1_B2M")),]

df_v2_2 <- df2

interacting_pair.x.new <- c()
j <- 0
for(x in as.vector(df_v2_2$interacting_pair.x)){
  j <- j + 1
  y <- paste0(str_split(x,"_")[[1]][2],"_",str_split(x,"_")[[1]][1])
  interacting_pair.x.new[j] <- y
}

CC.x.new <- c()
j <- 0
for(x in as.vector(df_v2_2$CC.x)){
  j <- j + 1
  y <- paste0(str_split(x,"\\.")[[1]][2],".",str_split(x,"\\.")[[1]][1])
  CC.x.new[j] <- y
}

df_v2_2$interacting_pair.x <- interacting_pair.x.new
df_v2_2$CC.x <- CC.x.new


df <- rbind(df1_use,df_v2_2)

df$CC.x <- paste0(df$group,"_",df$CC.x)

xlab <- c()
i <- 0
for(x in as.vector(df$CC.x)){
  i <- i + 1
  y <- paste0("<span style='color:blue'>",gsub("Glioma_|PBMC_","",str_split(x,"\\.")[[1]][1]),"</span><span style='color:black'> vs </span></span><span style='color:red'>",str_split(x,"\\.")[[1]][2],"</span>")
  xlab[i] <- y
}


ylab <- c()
j <- 0
for(x in as.vector(df$interacting_pair.x)){
  j <- j + 1
  if(x %in% c("CD3G_B2M","CD3D_B2M","C3AR1_C3","CCR4_CCL3","CCR5_CCL3","CCR5_CCL3L1")){
    y <- paste0("<span style='color:blue'>",str_split(x,"_")[[1]][2],"</span><span style='color:black'>-</span></span><span style='color:red'>",str_split(x,"_")[[1]][1],"</span>")
  }else{
    y <- paste0("<span style='color:blue'>",str_split(x,"_")[[1]][1],"</span><span style='color:black'>-</span></span><span style='color:red'>",str_split(x,"_")[[1]][2],"</span>")
  }
  ylab[j] <- y
}

df$xlab <- xlab
df$ylab <- ylab
df$log2means <- log2(df$means)

#df$ylab <- factor(df$ylab,levels=c(df$ylab[c(10,20,1,9,43)]))
df$ylab <- factor(df$ylab,levels=c(df$ylab[c(43,9,1,20,10)]))

breaks <- c(-2,-1,0,1,2)
pdf("Treg_vs_Myeloid_LR_v1.pdf",width=9,height=3.3)
ggplot(df, aes(x = xlab, y = ylab))+geom_point(aes(size = -log10(pvals+0.0001),color = log2means)) +
geom_vline(aes_string(xintercept = "xlab"),lwd = 0.2,colour = "grey",linetype="dashed") +
geom_hline(aes_string(yintercept = "ylab"),lwd = 0.2,colour = "grey",linetype="dashed")+scale_colour_gradientn(colors = c("#000000","#0000FF","#FFFF00","#FF0000"),breaks = breaks, labels = format(breaks),limits=c(-2,2)) +facet_grid(facets =  ~ group, space = "free", scales = "free")+theme_classic()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1))+theme(axis.text.y = element_markdown(), axis.text.x = element_markdown(angle = 90,hjust = 1,vjust = 0.5)) +
guides(size = guide_legend(title = "-log10(pvalue+1e-04)"),
       color = guide_colorbar(title = "Log2 mean (Molecule1, Molecule2)"))
dev.off()



#Others vs Treg interaction
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/cellphonedb/V1/compare_GBM_vs_PBMC")

#selection
GBMmeans <- read.csv("GBM_means.txt",header=T,sep="\t")
GBMpvals <- read.csv("GBM_pvalues.txt",header=T,sep="\t")

GBMmeans %>% dplyr::select("interacting_pair",starts_with(c("Treg")))  %>% reshape2::melt() -> GBMmeansdf

colnames(GBMmeansdf)<- c("interacting_pair","CC","means")

GBMpvals %>% dplyr::select("interacting_pair",starts_with(c("Treg")))%>% reshape2::melt()-> GBMpvalsdf

colnames(GBMpvalsdf)<- c("interacting_pair","CC","pvals")

GBMpvalsdf$joinlab<- paste0(GBMpvalsdf$interacting_pair,"_",GBMpvalsdf$CC)
GBMmeansdf$joinlab<- paste0(GBMmeansdf$interacting_pair,"_",GBMmeansdf$CC)
GBMpldf <- merge(GBMpvalsdf,GBMmeansdf,by = "joinlab")

GBMpldf$CC.x <- gsub("MDM\\.","MDM-",GBMpldf$CC.x)
GBMpldf$CC.x <- gsub("MG\\.","MG-",GBMpldf$CC.x)

GBMpldf2 <- GBMpldf[-grep("B$|NK$|Treg$",GBMpldf$CC.x),]

GBMpldf2 %>% filter(pvals < 1) %>% filter(means >1) -> df_v1

#df_v1$CC.x <- paste0("GBM_",df_v1$CC.x)
df_v1$group <- "Glioma"

###################
PBMCmeans <- read.csv("PBMC_means.txt",header=T,sep="\t")
PBMCpvals <- read.csv("PBMC_pvalues.txt",header=T,sep="\t")

PBMCmeans %>% dplyr::select("interacting_pair",starts_with(c("Treg")))  %>% reshape2::melt() -> PBMCmeansdf

colnames(PBMCmeansdf)<- c("interacting_pair","CC","means")

PBMCpvals %>% dplyr::select("interacting_pair",starts_with(c("Treg"))) %>% reshape2::melt()-> PBMCpvalsdf

colnames(PBMCpvalsdf)<- c("interacting_pair","CC","pvals")

PBMCpvalsdf$joinlab<- paste0(PBMCpvalsdf$interacting_pair,"_",PBMCpvalsdf$CC)
PBMCmeansdf$joinlab<- paste0(PBMCmeansdf$interacting_pair,"_",PBMCmeansdf$CC)
PBMCpldf <- merge(PBMCpvalsdf,PBMCmeansdf,by = "joinlab")

PBMCpldf2 <- PBMCpldf[-grep("B$|NK$|Treg$",PBMCpldf$CC.x),]

PBMCpldf2 %>% filter(pvals < 1) %>% filter(means >1) -> df_v2

#df_v1$CC.x <- paste0("GBM_",df_v1$CC.x)
df_v2$group <- "PBMC"

df1 <- rbind(df_v1,df_v2)

df1_use <- df1[which(df1$interacting_pair.x %in% c("CD3G_B2M","CD3D_B2M","CTLA4_CD86")),]


df_v1_2 <- df1_use

interacting_pair.x.new <- c()
j <- 0
for(x in as.vector(df_v1_2$interacting_pair.x)){
  j <- j + 1
  y <- paste0(str_split(x,"_")[[1]][2],"_",str_split(x,"_")[[1]][1])
  interacting_pair.x.new[j] <- y
}

CC.x.new <- c()
j <- 0
for(x in as.vector(df_v1_2$CC.x)){
  j <- j + 1
  y <- paste0(str_split(x,"\\.")[[1]][2],".",str_split(x,"\\.")[[1]][1])
  CC.x.new[j] <- y
}

df_v1_2$interacting_pair.x <- interacting_pair.x.new
df_v1_2$CC.x <- CC.x.new


GBMmeans <- read.csv("GBM_means.txt",header=T,sep="\t")
GBMpvals <- read.csv("GBM_pvalues.txt",header=T,sep="\t")

GBMmeans %>% dplyr::select("interacting_pair",starts_with(c("Monocyte","CD4","CD8","MDM","MG","DC")))  %>% reshape2::melt() -> GBMmeansdf

colnames(GBMmeansdf)<- c("interacting_pair","CC","means")

GBMpvals %>% dplyr::select("interacting_pair",starts_with(c("Monocyte","CD4","CD8","MDM","MG","DC")))%>% reshape2::melt()-> GBMpvalsdf

colnames(GBMpvalsdf)<- c("interacting_pair","CC","pvals")

GBMpvalsdf$joinlab<- paste0(GBMpvalsdf$interacting_pair,"_",GBMpvalsdf$CC)
GBMmeansdf$joinlab<- paste0(GBMmeansdf$interacting_pair,"_",GBMmeansdf$CC)
GBMpldf <- merge(GBMpvalsdf,GBMmeansdf,by = "joinlab")

GBMpldf$CC.x <- gsub("MDM\\.","MDM-",GBMpldf$CC.x)
GBMpldf$CC.x <- gsub("MG\\.","MG-",GBMpldf$CC.x)

GBMpldf2 <- GBMpldf[grep("Treg$",GBMpldf$CC.x),]

GBMpldf2 %>% filter(pvals < 1) %>% filter(means >1) -> df_v1

#df_v1$CC.x <- paste0("GBM_",df_v1$CC.x)
df_v1$group <- "Glioma"

###################
PBMCmeans <- read.csv("PBMC_means.txt",header=T,sep="\t")
PBMCpvals <- read.csv("PBMC_pvalues.txt",header=T,sep="\t")

PBMCmeans %>% dplyr::select("interacting_pair",starts_with(c("Monocyte","CD4","CD8","MDM","DC")))  %>% reshape2::melt() -> PBMCmeansdf

colnames(PBMCmeansdf)<- c("interacting_pair","CC","means")

PBMCpvals %>% dplyr::select("interacting_pair",starts_with(c("Monocyte","CD4","CD8","MDM","DC"))) %>% reshape2::melt()-> PBMCpvalsdf

colnames(PBMCpvalsdf)<- c("interacting_pair","CC","pvals")

PBMCpvalsdf$joinlab<- paste0(PBMCpvalsdf$interacting_pair,"_",PBMCpvalsdf$CC)
PBMCmeansdf$joinlab<- paste0(PBMCmeansdf$interacting_pair,"_",PBMCmeansdf$CC)
PBMCpldf <- merge(PBMCpvalsdf,PBMCmeansdf,by = "joinlab")

PBMCpldf2 <- PBMCpldf[grep("Treg$",PBMCpldf$CC.x),]

PBMCpldf2 %>% filter(pvals < 1) %>% filter(means >1) -> df_v2

#df_v1$CC.x <- paste0("GBM_",df_v1$CC.x)
df_v2$group <- "PBMC"

df1 <- rbind(df_v1,df_v2)

df2 <- df1[which(df1$interacting_pair.x %in% c("CCL3_CCR5","CCL3_CCR4","CCL3L1_CCR5")),]

df <- rbind(df_v1_2,df2)

df$CC.x <- paste0(df$group,"_",df$CC.x)

xlab <- c()
i <- 0
for(x in as.vector(df$CC.x)){
  i <- i + 1
  y <- paste0("<span style='color:blue'>",gsub("Glioma_|PBMC_","",str_split(x,"\\.")[[1]][1]),"</span><span style='color:black'> vs </span></span><span style='color:red'>",str_split(x,"\\.")[[1]][2],"</span>")
  xlab[i] <- y
}


ylab <- c()
j <- 0
for(x in as.vector(df$interacting_pair.x)){
  j <- j + 1
  if(x %in% c("CD3G_B2M","CD3D_B2M","C3AR1_C3","CCR4_CCL3","CCR5_CCL3","CCR5_CCL3L1")){
    y <- paste0("<span style='color:red'>",str_split(x,"_")[[1]][1],"</span><span style='color:black'>-</span></span><span style='color:blue'>",str_split(x,"_")[[1]][2],"</span>")
  }else{
    y <- paste0("<span style='color:blue'>",str_split(x,"_")[[1]][1],"</span><span style='color:black'>-</span></span><span style='color:red'>",str_split(x,"_")[[1]][2],"</span>")
  }
  ylab[j] <- y
}

df$xlab <- xlab
df$ylab <- ylab
df$log2means <- log2(df$means)

#df$ylab <- factor(df$ylab,levels=c(df$ylab[c(1,15,57,29,65,69)]))
df$ylab <- factor(df$ylab,levels=c(df$ylab[c(69,65,29,57,15,1)]))

breaks <- c(-2,-1,0,1,2)
pdf("Myeloid_vs_Treg_LR_v2.pdf",width=9,height=3.5)
ggplot(df, aes(x = xlab, y = ylab))+geom_point(aes(size = -log10(pvals+0.0001),color = log2means)) +
geom_vline(aes_string(xintercept = "xlab"),lwd = 0.2,colour = "grey",linetype="dashed") +
geom_hline(aes_string(yintercept = "ylab"),lwd = 0.2,colour = "grey",linetype="dashed")+scale_colour_gradientn(colors = c("#000000","#0000FF","#FFFF00","#FF0000"),breaks = breaks, labels = format(breaks),limits=c(-2,2)) +facet_grid(facets =  ~ group, space = "free", scales = "free")+theme_classic()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1))+theme(axis.text.y = element_markdown(), axis.text.x = element_markdown(angle = 90,hjust = 1,vjust = 0.5)) +
guides(size = guide_legend(title = "-log10(pvalue+1e-04)"),
       color = guide_colorbar(title = "Log2 mean (Molecule1, Molecule2)"))
dev.off()


#Glioma MDM vs MG
setwd("/data/activate_data/ranxiaojuan/ranxiaojuan_GBM_vs_PBMC_analysis/9sample_v2/cellphonedb/V1/compare_GBM_vs_PBMC")

GBMmeans <- read.csv("GBM_means.txt",header=T,sep="\t")
GBMpvals <- read.csv("GBM_pvalues.txt",header=T,sep="\t")

GBMmeans %>% dplyr::select("interacting_pair",starts_with(c("MDM.1","MDM.2","MDM.3")))  %>% reshape2::melt() -> GBMmeansdf

colnames(GBMmeansdf)<- c("interacting_pair","CC","means")

GBMpvals %>% dplyr::select("interacting_pair",starts_with(c("MDM.1","MDM.2","MDM.3")))%>% reshape2::melt()-> GBMpvalsdf

colnames(GBMpvalsdf)<- c("interacting_pair","CC","pvals")

GBMpvalsdf$joinlab<- paste0(GBMpvalsdf$interacting_pair,"_",GBMpvalsdf$CC)
GBMmeansdf$joinlab<- paste0(GBMmeansdf$interacting_pair,"_",GBMmeansdf$CC)
GBMpldf <- merge(GBMpvalsdf,GBMmeansdf,by = "joinlab")


GBMpldf2 <- GBMpldf[grep("MG.1$|MG.2$|MG.3",GBMpldf$CC.x),]
GBMpldf2$CC.x <- gsub("MDM\\.","MDM-",GBMpldf2$CC.x)
GBMpldf2$CC.x <- gsub("MG\\.","MG-",GBMpldf2$CC.x)

GBMpldf2 %>% filter(pvals < 1) %>% filter(means >1) -> df_v1

df1 <- df_v1[which(df_v1$interacting_pair.x %in% c("LDLR_APOE","PTPRC_LGALS1","APOE_LILRB4","APOE_TREM2","MIF_CD44","CCL4_CCR5","CCL4_CCR1","SPP1_CD44","SPP1_ITGA4")),]


interacting_pair.x.new <- c()
CC.x.new <- c()
j <- 0
for(x in as.vector(df1$interacting_pair.x)){
  j <- j + 1
  if(x %in% c("LDLR_APOE","PTPRC_LGALS1")){
    a <- paste0(str_split(x,"_")[[1]][2],"_",str_split(x,"_")[[1]][1])
  }else{
    a <- x
  }
  interacting_pair.x.new[j] <- a
}

i <- 0
for(x in as.vector(df1$interacting_pair.x)){
  i <- i + 1
  if (x %in% c("LDLR_APOE","PTPRC_LGALS1")){
      b <- paste0(str_split(df1$CC.x[i],"\\.")[[1]][2],".",str_split(df1$CC.x[i],"\\.")[[1]][1])
    }else{
      b <- df1$CC.x[i]
    }
  CC.x.new[i] <- b
}


df1$CC.x <- CC.x.new
df1$interacting_pair.x <- interacting_pair.x.new


GBMmeans <- read.csv("GBM_means.txt",header=T,sep="\t")
GBMpvals <- read.csv("GBM_pvalues.txt",header=T,sep="\t")

GBMmeans %>% dplyr::select("interacting_pair",starts_with(c("MG.1","MG.2","MG.3")))  %>% reshape2::melt() -> GBMmeansdf

colnames(GBMmeansdf)<- c("interacting_pair","CC","means")

GBMpvals %>% dplyr::select("interacting_pair",starts_with(c("MG.1","MG.2","MG.3")))%>% reshape2::melt()-> GBMpvalsdf

colnames(GBMpvalsdf)<- c("interacting_pair","CC","pvals")

GBMpvalsdf$joinlab<- paste0(GBMpvalsdf$interacting_pair,"_",GBMpvalsdf$CC)
GBMmeansdf$joinlab<- paste0(GBMmeansdf$interacting_pair,"_",GBMmeansdf$CC)
GBMpldf <- merge(GBMpvalsdf,GBMmeansdf,by = "joinlab")


GBMpldf2 <- GBMpldf[grep("MDM.1$|MDM.2$|MDM.3",GBMpldf$CC.x),]
GBMpldf2$CC.x <- gsub("MDM\\.","MDM-",GBMpldf2$CC.x)
GBMpldf2$CC.x <- gsub("MG\\.","MG-",GBMpldf2$CC.x)

GBMpldf2 %>% filter(pvals < 1) %>% filter(means >1) -> df_v1

df2 <- df_v1[which(df_v1$interacting_pair.x %in% c("APOE_TREM2","APOE_LILRB4","C1QA_CR1","VIM_CD44","IL1B_SIGIRR","HMGB1_CXCR4","CCL4_CCR1","CCL4_CCR5","C5AR1_RPS19","CD74_MIF")),]


interacting_pair.x.new <- c()
CC.x.new <- c()
j <- 0
for(x in as.vector(df2$interacting_pair.x)){
  j <- j + 1
  if(x %in% c("C5AR1_RPS19","CD74_MIF")){
    a <- paste0(str_split(x,"_")[[1]][2],"_",str_split(x,"_")[[1]][1])
  }else{
    a <- x
  }
  interacting_pair.x.new[j] <- a
}

i <- 0
for(x in as.vector(df2$interacting_pair.x)){
  i <- i + 1
  if (x %in% c("C5AR1_RPS19","CD74_MIF")){
      b <- paste0(str_split(df2$CC.x[i],"\\.")[[1]][2],".",str_split(df2$CC.x[i],"\\.")[[1]][1])
    }else{
      b <- df2$CC.x[i]
    }
  CC.x.new[i] <- b
}


df2$CC.x <- CC.x.new
df2$interacting_pair.x <- interacting_pair.x.new

df_all <- rbind(df1,df2) 

df <- df_all[grep("^MDM",df_all$CC.x),]


xlab <- c()
i <- 0
for(x in as.vector(df$CC.x)){
  i <- i + 1
  y <- paste0("<span style='color:blue'>",str_split(x,"\\.")[[1]][1],"</span><span style='color:black'>-</span></span><span style='color:red'>",str_split(x,"\\.")[[1]][2],"</span>")
  xlab[i] <- y
}


ylab <- c()
j <- 0
for(x in as.vector(df$interacting_pair.x)){
  j <- j + 1
  y <- paste0("<span style='color:blue'>",str_split(x,"_")[[1]][1],"</span><span style='color:black'>-</span></span><span style='color:red'>",str_split(x,"_")[[1]][2],"</span>")
  ylab[j] <- y
}

df$xlab <- xlab
df$ylab <- ylab
df$log2means <- log2(df$means)
#df$ylab <- factor(df$ylab,levels=c(df$ylab[c(9,1,18,24,37,39,28,31,34)]))
df$ylab <- factor(df$ylab,levels=c(df$ylab[c(34,31,28,39,37,24,18,1,9)]))


breaks <- c(-2,-1,0,1,2)
pdf("MDM_vs_MG_LR.pdf",width=6,height=3)
ggplot(df, aes(x = xlab, y = ylab))+geom_point(aes(size = -log10(pvals+0.0001),color = log2means)) +
geom_vline(aes_string(xintercept = "xlab"),lwd = 0.2,colour = "grey",linetype="dashed") +
geom_hline(aes_string(yintercept = "ylab"),lwd = 0.2,colour = "grey",linetype="dashed")+scale_colour_gradientn(colors = c("#000000","#0000FF","#FFFF00","#FF0000"),breaks = breaks, labels = format(breaks),limits=c(-2,2))+theme_classic()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5))+theme(axis.text.y = element_markdown(), axis.text.x = element_markdown(angle = 90,hjust = 1,vjust = 0.5)) +
guides(size = guide_legend(title = "-log10(pvalue+1e-04)"),
       color = guide_colorbar(title = "Log2 mean (Molecule1, Molecule2)"))
dev.off()




df <- df_all[grep("^MG",df_all$CC.x),]


xlab <- c()
i <- 0
for(x in as.vector(df$CC.x)){
  i <- i + 1
  y <- paste0("<span style='color:blue'>",str_split(x,"\\.")[[1]][1],"</span><span style='color:black'>-</span></span><span style='color:red'>",str_split(x,"\\.")[[1]][2],"</span>")
  xlab[i] <- y
}


ylab <- c()
j <- 0
for(x in as.vector(df$interacting_pair.x)){
  j <- j + 1
  y <- paste0("<span style='color:blue'>",str_split(x,"_")[[1]][1],"</span><span style='color:black'>-</span></span><span style='color:red'>",str_split(x,"_")[[1]][2],"</span>")
  ylab[j] <- y
}

df$xlab <- xlab
df$ylab <- ylab
df$log2means <- log2(df$means)
#df$ylab <- factor(df$ylab,levels=c(df$ylab[c(10,19,1,31,44,7,37,41,28,50)]))
df$ylab <- factor(df$ylab,levels=c(df$ylab[c(50,28,41,37,7,44,31,1,19,10)]))


breaks <- c(-2,-1,0,1,2)
pdf("MG_vs_MDM_LR.pdf",width=6,height=3.3)
ggplot(df, aes(x = xlab, y = ylab))+geom_point(aes(size = -log10(pvals+0.0001),color = log2means)) +
geom_vline(aes_string(xintercept = "xlab"),lwd = 0.2,colour = "grey",linetype="dashed") +
geom_hline(aes_string(yintercept = "ylab"),lwd = 0.2,colour = "grey",linetype="dashed")+scale_colour_gradientn(colors = c("#000000","#0000FF","#FFFF00","#FF0000"),breaks = breaks, labels = format(breaks),limits=c(-2,2))+theme_classic()+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=0.5))+theme(axis.text.y = element_markdown(), axis.text.x = element_markdown(angle = 90,hjust = 1,vjust = 0.5)) +
guides(size = guide_legend(title = "-log10(pvalue+1e-04)"),
       color = guide_colorbar(title = "Log2 mean (Molecule1, Molecule2)"))
dev.off()
