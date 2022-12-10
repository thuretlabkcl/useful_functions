
library(genefilter)
library(tidyverse)

# for DDS object (DESeq output object)
plotPCA.any.dds <- function (object, pc.x, pc.y, labels = FALSE, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  require(genefilter)
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  
  ## Select the PCAs and percentVar that you like instead of 1 and 2
  d <- data.frame(PCX = pca$x[, pc.x], PCY = pca$x[, pc.y], group = group, 
                  intgroup.df, sample = rownames(colData(object)))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pc.x:pc.y]
    return(d)
  }
  if (labels) {
    g <- ggplot(data = d, aes_string(x = "PCX", y = "PCY", color = "group")) + 
      geom_point(size = 3) + 
      xlab(paste0("PC", pc.x, ": ", round(percentVar[pc.x] * 100), "% Variance")) + 
      ylab(paste0("PC", pc.y,  ": ", round(percentVar[pc.y] * 100), "% Variance")) + 
      coord_fixed() + 
      geom_text_repel(size=3, aes_string(label = "sample"), color = "black") +
      scale_color_discrete(name = intgroup)
    return (g)
  } else {
    p <- ggplot(data = d, aes_string(x = "PCX", y = "PCY", color = "group")) + 
      geom_point(size = 3) + 
      xlab(paste0("PC", pc.x, ": ", round(percentVar[pc.x] * 100), "% Variance")) + 
      ylab(paste0("PC", pc.y,  ": ", round(percentVar[pc.y] * 100), "% Variance")) + 
      coord_fixed() + 
      scale_color_discrete(name = intgroup)
    return(p)
  }
}

# or for ExpressionSet object 

plotPCA.any.eset <- function (eset, pc.x, pc.y, labels = FALSE, intgroup, ntop = 500, returnData = FALSE) 
{
  require(genefilter)
  rv <- rowVars(exprs(eset))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(exprs(eset)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(pData(eset)))) {
    stop("the argument 'intgroup' should specify columns of phenotype data in expression set object")
  }
  intgroup.df <- as.data.frame(pData(eset)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    pData(eset)[[intgroup]]
  }
  
  ## Select the PCAs and percentVar that you like instead of 1 and 2
  d <- data.frame(PCX = pca$x[, pc.x], PCY = pca$x[, pc.y], group = group, 
                  intgroup.df, sample = rownames(pData(eset)))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pc.x:pc.y]
    return(d)
  }
  
  ## CREATE PLOTS ###
  require(ggplot2)
  require(ggrepel)
  
  if (labels) {
    g <- ggplot(data = d, aes_string(x = "PCX", y = "PCY", color = "group")) + 
      geom_point(size = 3) + 
      xlab(paste0("PC", pc.x, ": ", round(percentVar[pc.x] * 100), "% Variance")) + 
      ylab(paste0("PC", pc.y,  ": ", round(percentVar[pc.y] * 100), "% Variance")) + 
      coord_fixed() + 
      geom_text_repel(size=3, aes_string(label = "sample"), color = "black") +
      scale_color_discrete(name = intgroup) +
      theme_bw() +
      theme(legend.text = element_text(size = 16), 
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14),
            legend.key.size = unit(1.5, "cm")) 
    return (g)
  } else {
    p <- ggplot(data = d, aes_string(x = "PCX", y = "PCY", color = "group")) + 
      geom_point(size = 3) + 
      xlab(paste0("PC", pc.x, ": ", round(percentVar[pc.x] * 100), "% Variance")) + 
      ylab(paste0("PC", pc.y,  ": ", round(percentVar[pc.y] * 100), "% Variance")) + 
      coord_fixed() + 
      scale_color_discrete(name = intgroup) +
      theme_bw() +
      theme(legend.text = element_text(size = 16), 
            legend.title = element_text(size = 18),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 14),
            legend.key.size = unit(1.5, "cm"))
    return(p)
  }
}
