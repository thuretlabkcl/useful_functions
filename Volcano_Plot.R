

### VOLCANO PLOT SEURAT ################################################################################

## READ ME ###

# minimum input is the dataframe output from FindMarkers() call in Seurat package. It needs a column called "gene"
# a column called "avg_log2FC" and a another column called "p_val_adj"


plot_volcano_Seurat <- function (results, n_up = 2, n_down = 2, bylogFC = FALSE, 
                          genelabels = TRUE, point.padding = 0.22, textsize = 2.8) 
{
  
  require(ggplot2)
  require(ggforce)
  require(tidyverse)
  
  res_1 <- data.frame(results) %>% 
    mutate(threshold = p_val_adj < 0.05) %>%
    arrange(p_val_adj)
  
  
  res_2 <- res_1 %>% # filter out rows containing no multiple corrected p value
    filter(!is.na(p_val_adj))
  
  
  res_col <- res_2 %>% # label genes by whether there are significantly up or down regulated
    mutate(up_down = ifelse(avg_log2FC < 0 & p_val_adj < 0.05, "Down", 
                            ifelse(avg_log2FC > 0 & p_val_adj < 0.05, "Up", "Not Sig"))) 
  
  
  res_3 <- res_col %>% # get df with only padj < 0.05
    filter(up_down == "Up" | up_down == "Down") 
  
  if (bylogFC) {
    res_4 <- res_3 %>% 
      arrange(desc(up_down), desc(avg_log2FC), p_val_adj) # get order of up genes by descending logFC
    
    res_4 <- res_4 %>%  # grab top n genes for labelling 
      mutate(genelabs = ifelse(res_4$gene %in% res_4$gene[1:n_up],
                               as.character(res_4$gene), "")) 
    res_5 <- res_4 %>% 
      arrange(up_down,avg_log2FC, p_val_adj) # get order of down genes by logFC lowest to highest
    
    res_6 <- res_5 %>% 
      mutate(genelabs = ifelse(res_5$gene %in% res_5$gene[1:n_down], 
                               as.character(res_5$gene), ""))
    
    res_labs <- full_join(res_4, res_6) # merge to get final genelabs column containing up and down gene labels
    
  } else {
    res_4 <- res_3 %>% 
      arrange(desc(up_down), p_val_adj) # get order of up genes sorted by padj value smallest to largest
    
    res_4 <- res_4 %>% 
      mutate(genelabs = ifelse(res_4$gene %in% res_4$gene[1:n_up],
                               as.character(res_4$gene), "")) 
    res_5 <- res_4 %>% 
      arrange(up_down, p_val_adj) # get order of down genes sorted by padj value smallest to largest
    
    res_6 <- res_5 %>% 
      mutate(genelabs = ifelse(res_5$gene %in% res_5$gene[1:n_down], 
                               as.character(res_5$gene), ""))
    
    res_labs <- full_join(res_4, res_6) # merge to get final genelabs column containing up and down gene labels
  }
  
  # CREATE VOLCANO PLOT
  
  require(ggplot2)
  require(ggrepel)
  require(RColorBrewer)
  require(ggforce)
  
  if (isTRUE(genelabels)) {
    plot <- ggplot(res_col, aes(x = avg_log2FC, y = -log10(p_val_adj), 
                                col = up_down)) + 
      geom_point(show.legend = F, alpha = 0.7) +
      theme_bw() +
      scale_color_brewer(palette = 2, type = "qual") +
      geom_hline(yintercept = -log10(0.05), color = "black", linetype = 2, size = 0.3) +
      geom_vline(xintercept = -0.32, color = "black", linetype = 2, size = 0.3) +
      geom_vline(xintercept = 0.32, color = "black", linetype = 2, size = 0.3) +
      geom_text_repel(data = res_labs, aes(label = genelabs), size = textsize, color = "black", point.padding = point.padding) +
      theme(axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 13)) +
      labs(x = "Log2 Fold Change", y = "-log10 (Adjusted P Value)")
    return(plot)
  } else {
    plot <- ggplot(res_col, aes(x = avg_log2FC, y = -log10(p_val_adj), 
                                col = up_down)) + 
      geom_point(show.legend = F, alpha = 0.7) +
      theme_bw() +
      scale_color_brewer(palette = 2, type = "qual") +
      geom_hline(yintercept = -log10(0.05), color = "black", linetype = 2, size = 0.3) +
      geom_vline(xintercept = -0.32, color = "black", linetype = 2, size = 0.3) +
      geom_vline(xintercept = 0.32, color = "black", linetype = 2, size = 0.3) +
      theme(axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 13)) +
      labs(x = "Log2 Fold Change", y = "-log10 (Adjusted P Value)")
    return(plot)
  }
}  




# n_up = number of upregulated genes to label
# n_down = number of downregulated genes to label
# bylogFC = whether to label the points according to log 2 Fold Change (TRUE), or padj significance level (FALSE)
# genelabels = turn text labelling of gene points on or off
# point.padding = degree of space between points and gene labels
# textsize = fontsize of gene labels




# e.g 
#volcano_plot(results = prog_res_shrunk, n_up = 8, n_down = 3, 
             #genelabels = T, bylogFC = T, point.padding = 0.3, textsize = 3)





### VOLCANO PLOT DESEQ2 ############################################################################################

 ## READ ME ###

# minimum input is a results object from a results() call on a dds object run through DESEQ()
# better input is a logFC shrunk results object to give better DE estimation


volcano_plot_DESeq2 <- function (results, n_up = 2, n_down = 2, bylogFC = FALSE, 
                          genelabels = TRUE, point.padding = 0.22, textsize = 2.8) 
  {
  
  require(dplyr)
  genes <- rownames(results)
  res_1 <- data.frame(results) %>% 
    mutate(threshold = padj < 0.05, 
           gene = genes) %>% # tibble deletes rownames so must create column to store
    arrange(padj)
  
  
  res_2 <- res_1 %>% # filter out rows containing no multiple corrected p value
    filter(!is.na(padj))
  
  
  res_col <- res_2 %>% # label genes by whether there are significantly up or down regulated
    mutate(up_down = ifelse(log2FoldChange < 0 & padj < 0.05, "Down", 
                            ifelse(log2FoldChange > 0 & padj < 0.05, "Up", "Not Sig"))) 
  
  
  res_3 <- res_col %>% # get df with only padj < 0.05
    filter(up_down == "Up" | up_down == "Down") 
  
  if (bylogFC) {
    res_4 <- res_3 %>% 
      arrange(desc(up_down), desc(log2FoldChange), padj) # get order of up genes by descending logFC
    
    res_4 <- res_4 %>%  # grab top n genes for labelling 
      mutate(genelabs = ifelse(res_4$gene %in% res_4$gene[1:n_up],
                               as.character(res_4$gene), "")) 
    res_5 <- res_4 %>% 
      arrange(up_down,log2FoldChange, padj) # get order of down genes by logFC lowest to highest
    
    res_6 <- res_5 %>% 
      mutate(genelabs = ifelse(res_5$gene %in% res_5$gene[1:n_down], 
                               as.character(res_5$gene), ""))
    
    res_labs <- full_join(res_4, res_6) # merge to get final genelabs column containing up and down gene labels
    
  } else {
    res_4 <- res_3 %>% 
      arrange(desc(up_down), padj) # get order of up genes sorted by padj value smallest to largest
      
    res_4 <- res_4 %>% 
      mutate(genelabs = ifelse(res_4$gene %in% res_4$gene[1:n_up],
                               as.character(res_4$gene), "")) 
    res_5 <- res_4 %>% 
      arrange(up_down, padj) # get order of down genes sorted by padj value smallest to largest
    
    res_6 <- res_5 %>% 
      mutate(genelabs = ifelse(res_5$gene %in% res_5$gene[1:n_down], 
                               as.character(res_5$gene), ""))
    
    res_labs <- full_join(res_4, res_6) # merge to get final genelabs column containing up and down gene labels
  }
 
  # CREATE VOLCANO PLOT
  
  require(ggplot2)
  require(ggrepel)
  require(RColorBrewer)
  require(ggforce)
  
  if (isTRUE(genelabels)) {
    plot <- ggplot(res_col, aes(x = log2FoldChange, y = -log10(padj), 
                                col = up_down)) + 
      geom_point(show.legend = F, alpha = 0.7) +
      theme_bw() +
      scale_color_brewer(palette = 2, type = "qual") +
      geom_hline(yintercept = -log10(0.05), color = "black", linetype = 2, size = 0.3) +
      geom_vline(xintercept = -0.32, color = "black", linetype = 2, size = 0.3) +
      geom_vline(xintercept = 0.32, color = "black", linetype = 2, size = 0.3) +
      geom_text_repel(data = res_labs, aes(label = genelabs), size = textsize, color = "black", point.padding = point.padding) +
      theme(axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.text = element_text(size = 13)) +
      labs(x = "Log2 Fold Change", y = "-log10 (Adjusted P Value)")
    return(plot)
    } else {
      plot <- ggplot(res_col, aes(x = log2FoldChange, y = -log10(padj), 
                                  col = up_down)) + 
        geom_point(show.legend = F, alpha = 0.7) +
        theme_bw() +
        scale_color_brewer(palette = 2, type = "qual") +
        geom_hline(yintercept = -log10(0.05), color = "black", linetype = 2, size = 0.3) +
        geom_vline(xintercept = -0.32, color = "black", linetype = 2, size = 0.3) +
        geom_vline(xintercept = 0.32, color = "black", linetype = 2, size = 0.3) +
        theme(axis.title.x = element_text(size = 16),
              axis.title.y = element_text(size = 16),
              axis.text = element_text(size = 13)) +
        labs(x = "Log2 Fold Change", y = "-log10 (Adjusted P Value)")
      return(plot)
    }
  }  



# n_up = number of upregulated genes to label
# n_down = number of downregulated genes to label
# bylogFC = whether to label the points according to log 2 Fold Change (TRUE), or padj significance level (FALSE)
# genelabels = turn text labelling of gene points on or off
# point.padding = degree of space between points and gene labels
# textsize = fontsize of gene labels




# e.g 
#volcano_plot(results = prog_res_shrunk, n_up = 8, n_down = 3, 
             # genelabels = T, bylogFC = T, point.padding = 0.3, textsize = 3)



### VOLCANO PLOT LIMMA ############################################################################################

  # This function has been adapted to handle a results dataframe from the output of a LIMMA analysis
  
  volcano_plot_limma <- function (results, n_up = 2, n_down = 2, bylogFC = FALSE, 
                                  genelabels = TRUE, point.padding = 0.22, textsize = 2.8) 
  {
    
    require(dplyr)
    genes <- results$GENE_SYMBOL
    res_1 <- data.frame(prolif_results) %>%
      rename(padj = adj.P.Val,log2FoldChange = logFC,p = P.Value) %>% 
      mutate(threshold = p < 0.05, 
             gene = genes) %>% # tibble deletes rownames so must create column to store
      arrange(p)
    
    
    res_2 <- res_1 %>% # filter out rows containing no multiple corrected p value
      filter(!is.na(p))
    
    
    res_col <- res_2 %>% # label genes by whether there are significantly up or down regulated
      mutate(up_down = ifelse(log2FoldChange < 0 & p < 0.05, "Down", 
                              ifelse(log2FoldChange > 0 & p < 0.05, "Up", "Not Sig"))) 
    
    
    res_3 <- res_col %>% # get df with only p < 0.05
      filter(up_down == "Up" | up_down == "Down") 
    
    if (bylogFC) {
      res_4 <- res_3 %>% 
        arrange(desc(up_down), desc(log2FoldChange), p) # get order of up genes by descending logFC
      
      res_4 <- res_4 %>%  # grab top n genes for labelling 
        mutate(genelabs = ifelse(res_4$gene %in% res_4$gene[1:n_up],
                                 as.character(res_4$gene), "")) 
      res_5 <- res_4 %>% 
        arrange(up_down,log2FoldChange, p) # get order of down genes by logFC lowest to highest
      
      res_6 <- res_5 %>% 
        mutate(genelabs = ifelse(res_5$gene %in% res_5$gene[1:n_down], 
                                 as.character(res_5$gene), ""))
      
      res_labs <- full_join(res_4, res_6) # merge to get final genelabs column containing up and down gene labels
      
    } else {
      res_4 <- res_3 %>% 
        arrange(desc(up_down), p) # get order of up genes sorted by p value smallest to largest
      
      res_4 <- res_4 %>% 
        mutate(genelabs = ifelse(res_4$gene %in% res_4$gene[1:n_up],
                                 as.character(res_4$gene), "")) 
      res_5 <- res_4 %>% 
        arrange(up_down, p) # get order of down genes sorted by p value smallest to largest
      
      res_6 <- res_5 %>% 
        mutate(genelabs = ifelse(res_5$gene %in% res_5$gene[1:n_down], 
                                 as.character(res_5$gene), ""))
      
      res_labs <- full_join(res_4, res_6) # merge to get final genelabs column containing up and down gene labels
    }
    
    # CREATE VOLCANO PLOT
    
    require(ggplot2)
    require(ggrepel)
    require(RColorBrewer)
    require(ggforce)
    
    if (isTRUE(genelabels)) {
      plot <- ggplot(res_col, aes(x = log2FoldChange, y = -log10(p), 
                                  col = up_down)) + 
        geom_point(show.legend = F, alpha = 0.7) +
        theme_bw() +
        scale_color_brewer(palette = 2, type = "qual") +
        geom_hline(yintercept = -log10(0.05), color = "black", linetype = 2, size = 0.3) +
        geom_vline(xintercept = -1, color = "black", linetype = 2, size = 0.3) +
        geom_vline(xintercept = 1, color = "black", linetype = 2, size = 0.3) +
        geom_text_repel(data = res_labs, aes(label = genelabs), size = textsize, color = "black", point.padding = point.padding) +
        theme(axis.title.x = element_text(size = 16),
              axis.title.y = element_text(size = 16),
              axis.text = element_text(size = 16)) +
        labs(x = "Log2 Fold Change", y = "-log10 (P Value)")
      return(plot)
    } else {
      plot <- ggplot(res_col, aes(x = log2FoldChange, y = -log10(p), 
                                  col = up_down)) + 
        geom_point(show.legend = F, alpha = 0.7) +
        theme_bw() +
        scale_color_brewer(palette = 2, type = "qual") +
        geom_hline(yintercept = -log10(0.05), color = "black", linetype = 2, size = 0.3) +
        geom_vline(xintercept = -1, color = "black", linetype = 2, size = 0.3) +
        geom_vline(xintercept = 1, color = "black", linetype = 2, size = 0.3) +
        theme(axis.title.x = element_text(size = 16),
              axis.title.y = element_text(size = 16),
              axis.text = element_text(size = 16)) +
        labs(x = "Log2 Fold Change", y = "-log10 (P Value)")
      return(plot)
    }
  }  

