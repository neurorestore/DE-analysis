recode_colnames = function(DE) {
  colnames(DE) %<>%
    fct_recode('p_val' = 'p.value',  ## DESeq2
               'p_val' = 'pvalue',  ## DESeq2
               'p_val' = 'p.value',  ## t/wilcox
               'p_val' = 'P.Value',  ## limma
               'p_val' = 'PValue'  , ## edgeR
               'p_val_adj' = 'padj', ## DESeq2/t/wilcox
               'p_val_adj' = 'adj.P.Val',      ## limma
               'p_val_adj' = 'FDR',            ## edgeER
               'avg_logFC' = 'log2FoldChange', ## DESEeq2
               'avg_logFC' = 'logFC', ## limma/edgeR
               'test_statistic' = 'stat', ## DESeq2
               'test_statistic' = 'F', ## edgeR
               'test_statistic' = 't', ## limma
               'test_statistic' = 'LR', ## edgeR LRT
               'test_statistic' = 'statistic' ## t
    ) %>%
    as.character()  
  return(DE)
}
