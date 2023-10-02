suppressPackageStartupMessages({
library("ComplexHeatmap")
library("circlize")
library("gridBase")
library("grid")
library("RColorBrewer")
library("zoo")
})

# 1. adjust on sscClust::ssc.plot.headmap
heatmap_sm = function (obj, assay.name = "exprs", out.prefix = NULL, ncell.downsample = NULL, 
  ave.by = NULL, columns = NULL, columns.order = NULL, gene.desc = NULL, 
  colSet = list(), pdf.width = 16, pdf.height = 15, do.scale = TRUE, 
  z.lo = -2.5, z.hi = 2.5, z.step = 1, exp.title = "Exp", 
  do.clustering.row = T, do.clustering.col = T, dend.col = FALSE, 
  dend.row = FALSE, clustering.distance = "spearman", clustering.method = "complete", 
  k.row = 1, k.col = 1, palette.name = NULL, annotation_legend_param = list(), 
  ann.bar.height = 1.5, mytitle = "", ...) 
{
  ## remove cell contain NA values
  flag <- assay(obj, assay.name) %>% t %>% complete.cases
  obj <- obj[,flag]
  ##
  if (!is.null(gene.desc) && ("Group" %in% colnames(gene.desc)) && 
    ("geneID" %in% colnames(gene.desc))) {
    obj <- obj[gene.desc$geneID, ]
  }
  if (!is.null(ncell.downsample) && ncell.downsample < ncol(obj)) {
    obj <- obj[, sample(seq_len(ncol(obj)), ncell.downsample)]
  }
  n <- nrow(obj)
  m <- ncol(obj)
  if (n < 3) {
    loginfo(sprintf("Too few genes: n=%s", n))
    return(NULL)
  }
  if (m < 3) {
    loginfo(sprintf("Too few samples: m=%s", m))
    return(NULL)
  }
  if (is.null(ave.by)) {
    obj <- ssc.assay.hclust(obj, assay.name = assay.name, 
      order.col = if (is.logical(dend.col) && FALSE == 
        dend.col) 
        do.clustering.col
      else FALSE, order.row = if (is.logical(dend.row) && 
        FALSE == dend.row) 
        do.clustering.row
      else FALSE, clustering.distance = "spearman", clustering.method = "complete", 
      k.row = 1, k.col = 1)
  }
  else {
    obj <- ssc.average.cell(obj, assay.name = assay.name, 
      column = ave.by, ret.type = "sce")
    columns <- intersect(ave.by, columns)
    columns.order <- intersect(ave.by, columns.order)
  }
  ha.col <- NULL
  annDF <- data.frame()
  if (!is.null(columns)) {
    if (!is.null(columns.order)) {
      obj <- ssc.order(obj, columns.order = columns.order)
    }
    annDF <- as.data.frame(colData(obj)[columns])
    if (length(colSet) == 0) {
      for (i in seq_along(columns)) {
        x <- columns[i]
        if (class(colData(obj)[, x]) == "numeric") {
          if (all(colData(obj)[, x] <= 1) && all(colData(obj)[, 
            x] >= 0)) {
            Y.level <- c(0, 1)
          }
          else {
            Y.level <- pretty(colData(obj)[, x], n = 8)
          }
          colSet[[x]] <- colorRamp2(seq(Y.level[1], 
            Y.level[length(Y.level)], length = 7), rev(brewer.pal(n = 7, 
            name = "RdYlBu")), space = "LAB")
          annotation_legend_param[[x]] <- list(color_bar = "continuous", 
            legend_direction = "horizontal", legend_width = unit(4, 
              "cm"), legend_height = unit(2, "cm"))
        }
        else {
          group.value <- sort(unique(colData(obj)[, 
            x]))
          colSet[[x]] <- structure(auto.colSet(length(group.value), 
            name = "Accent"), names = group.value)
        }
      }
    }
    g.show.legend <- T
    ha.col <- ComplexHeatmap::HeatmapAnnotation(df = annDF, 
      col = colSet, show_legend = g.show.legend, simple_anno_size = unit(ann.bar.height, 
        "cm"), annotation_legend_param = annotation_legend_param)
  }
  obj <- ssc.order(obj, columns.order = NULL, gene.desc = gene.desc)
  dat.plot <- as.matrix(assay(obj, assay.name))
  rownames(dat.plot) <- unname(rowData(obj)$display.name)
  if (do.scale) {
    rowM <- rowMeans(dat.plot, na.rm = T)
    rowSD <- apply(dat.plot, 1, sd, na.rm = T)
    dat.plot <- sweep(dat.plot, 1, rowM)
    dat.plot <- sweep(dat.plot, 1, rowSD, "/")
    if (!is.null(z.lo)) {
      dat.plot[dat.plot < z.lo] <- z.lo
    }
    if (!is.null(z.hi)) {
      dat.plot[dat.plot > z.hi] <- z.hi
    }
  }
  else {
    tmp.var <- pretty((dat.plot), n = 8)
    if (is.null(z.lo)) {
      z.lo <- tmp.var[1]
    }
    if (is.null(z.hi)) {
      z.hi <- tmp.var[length(tmp.var)]
    }
    if (is.null(z.step)) {
      z.step <- tmp.var[2] - tmp.var[1]
    }
  }
  if (!is.null(out.prefix)) {
    pdf(sprintf("%s.pdf", out.prefix), width = pdf.width, 
      height = pdf.height)
  }
  par(mar = c(4, 12, 4, 4))
  plot.new()
  title(main = mytitle, cex.main = 2)
  vps <- baseViewports()
  pushViewport(vps$inner, vps$figure, vps$plot)
  if (is.null(palette.name)) {
    exp.palette <- rev(brewer.pal(n = 7, name = ifelse(do.scale, 
      "RdBu", "RdYlBu")))
  }
  else {
	if (palette.name=="blue2green2red"){
		bks <- seq(-3.1, 3.1, by=0.1)
		exp.palette <- blue2green2red(length(bks) - 1)
	}
	else{
		exp.palette <- rev(brewer.pal(n = 7, name = palette.name))
	}
  }
  # smooth value
  dat.plot = t(apply(dat.plot, 1, function(x){rollmean(x, 50, fill="extend")}))
  #
  ht <- ComplexHeatmap::Heatmap(dat.plot, name = exp.title, 
    col = colorRamp2(seq(z.lo, z.hi, length = 100), colorRampPalette(exp.palette)(100), 
      space = "LAB"), column_dend_height = unit(6, "cm"), 
    row_dend_width = unit(1, "cm"), column_names_gp = grid::gpar(fontsize = 12 * 
      28/max(m, 32)), row_names_gp = grid::gpar(fontsize = 18 * 
      28/max(n, 32)), show_heatmap_legend = T, row_names_max_width = unit(10, 
      "cm"), cluster_columns = dend.col, cluster_rows = dend.row, 
    row_dend_reorder = FALSE, column_dend_reorder = FALSE, 
    heatmap_legend_param = list(grid_width = unit(0.8, "cm"), 
      grid_height = unit(0.8, "cm"), at = seq(z.lo, z.hi, 
        z.step), title_gp = grid::gpar(fontsize = 14, 
        fontface = "bold"), label_gp = grid::gpar(fontsize = 12), 
      color_bar = "continuous"), top_annotation = ha.col, 
    ...)
  ComplexHeatmap::draw(ht, newpage = FALSE)
  if (!is.null(ha.col)) {
    for (i in seq_along(names(ha.col@anno_list))) {
      ComplexHeatmap::decorate_annotation(names(ha.col@anno_list)[i], 
        {
          grid.text(names(ha.col@anno_list)[i], unit(-4, 
            "mm"), gp = grid::gpar(fontsize = 14), just = "right")
        })
    }
  }
  if (!is.null(out.prefix)) {
    dev.off()
  }
  return(ht)
}



## 2. new functions
convertGeneID = function(.genes, .from="SYMBOL", .to="ENSEMBL"){
        get = select(org.Hs.eg.db, keys=.genes, keytype=.from, columns=.to)
        get = get[,2]
        get = get[!is.na(get)]
        return(get)
}

convertDPTperc = function(.sce){
  .sce = .sce[,order(.sce$dpt_pseudotime)]
  .sce$dpt_order = 1:ncol(.sce)
  .sce$dpt_order_perc = .sce$dpt_order / max(.sce$dpt_order)
  return(.sce)
}

testOneGene = function(.gene, .sce){
  dat = data.frame(exp=as.vector(assay(.sce[.gene,],"exprs")), time=.sce$dpt_order_perc)
  dat = dat[!is.na(dat$exp),]
  #fit = gam(exp ~ lo(time), data=dat)
  fit = gam(exp ~ time, data=dat)
  summ = summary(fit)
  coef = unname(coef(fit)[2])
  p = summ[4][[1]][1,5]
  m_sq = summ[4][[1]][1,3]
  df = data.frame(geneSymbol=.gene, mean_sq=m_sq, coef=coef, pval=p)
  return(df)
}

testProcess = function(.sce, .tfs){
  # correaltion test
  gam.test = ldply( rownames(.sce), testOneGene, .sce, .parallel=T)
  colnames(gam.test) = c("geneSymbol","mean_sq", "coef", "pval")
  gam.test$geneSymbol = as.character(gam.test$geneSymbol)
  rownames(gam.test) = gam.test$geneSymbol
  gam.test = gam.test[order(gam.test$pval,decreasing=F),]
  gam.test$adj.p = p.adjust(gam.test$pval, method="BH")
  gam.test$is.TF = ifelse(as.character(gam.test$geneSymbol) %in% .tfs, T, F)
  return(gam.test)
}

plotONEgene = function(.gene, .sce, .colSet){
  dat = data.frame(expression=as.vector(assay(.sce[.gene,],"exprs")), 
                   pseudotime=.sce$dpt_order_perc,
                   meta.cluster=.sce$meta.cluster)
  dat = dat[!is.na(dat$expression),]
  #
  model = loess(as.numeric(dat$expression) ~ as.numeric(dat$pseudotime))
  Y = predict(model)
  X = dat$pseudotime
  dY = diff(Y)/diff(X)
  dY = rollmean(dY, 30, fill="extend")  # smooth using 30 size window
  dX = rowMeans(embed(X,2))
  new.df = data.frame(dX=dX, dY=dY)
  #
  y.cutoff=1
  breaks = matrix(NA, nrow=0, ncol=2)
  last.start = NA
  last.end  = NA
  while(1){
    point.start = new.df[dY>=y.cutoff & (is.na(last.end) | dX>last.end), "dX"] %>% head(1)
    point.end = new.df[dY<y.cutoff & dX>point.start, "dX"] %>% head(1)
    # find a break
    if(length(point.start)>0){
      if(length(point.end)==0){
        # to the end
        point.end = max(dX)
      }
      breaks = rbind(breaks, c(point.start,point.end))
      last.start = point.start
      last.end = point.end
    }else{
        break
    }
  }
  breaks = as.data.frame(breaks)
  colnames(breaks) = c("X1","X2")
  #breaks = breaks[which.max(breaks$X2-breaks$X1),]  # only keep the biggest one
  #
  #adj.P = .gam[.gene, "adj.p"]
  #info = sprintf(sprintf("%s\nP = %.3f", .gene, adj.P))
  info = .gene
  #
  p = ggplot(dat, aes(x=pseudotime, y=expression))
  if (nrow(breaks)!=0){
    for (i in 1:nrow(breaks)){
      p = p + ggplot2::geom_rect(xmin=breaks[i,"X1"], xmax=breaks[i,"X2"], ymin=-Inf, ymax=Inf, fill="lightgrey")
    }
  }
  p = p + geom_point_rast(aes(color=meta.cluster),shape=16, size=1,  raster.dpi=300) +
          geom_smooth(aes(color=NULL), color="black", size=1, method="loess", formula=y~x, se=F) +
          geom_line(data=new.df,aes(x=dX,y=dY), lty=5, size=1) +
          geom_hline(yintercept=y.cutoff, color="red", lty="dashed", size=0.6) +
          ggtitle(info) +
          theme_classic2() +
          scale_color_manual(values=.colSet$meta.cluster[unique(dat$meta.cluster)]) +
          xlab("dpt_order_perc") + 
          coord_cartesian(ylim=c(-1,2.5))
  return(p)
}

hyper.test = function(set1, set2, bg){
    over = intersect(set1,set2)
    q = length(over) -1
    m = length(set2)
    n = length(bg) - length(set2)
    k = length(set1)
    p.val = phyper(q, m, n, k, lower.tail=F)
    return(p.val)
}
