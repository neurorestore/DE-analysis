## function to score the overlap between bulk and single-cell DE
calculate_overlap = function(bulk_de, sc_de,
                             method = c('fcc', 'aucc'),
                             k = NULL,
                             cor_method = c('pearson', 'spearman')) {
  method = match.arg(method)
  cor__method = match.arg(cor_method)
  if (method == 'fcc' & is.na(cor_method)) {
    stop("if using method='fcc', you must set cor_method (pearson/spearman)")
  }
  if (method == 'aucc' & is.na(k)) {
    stop("If using method='aucc', you must set k")
  }
  
  # double check column names
  colnames(bulk_de) %<>%
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
  colnames(sc_de) %<>%
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
               'test_statistic' = 'LR', ## edgeR LRT
               'test_statistic' = 't', ## limma
               'test_statistic' = 'statistic' ## t
    ) %>%
    as.character()
  
  # remove NAs
  sc_de %<>% filter(!is.na(p_val), !is.na(p_val_adj), !is.na(test_statistic))
  bulk_de %<>% filter(!is.na(p_val), !is.na(p_val_adj), !is.na(test_statistic))
  
  # replace p=0 with minimum p-value
  sc_de_min = min(sc_de$p_val_adj[sc_de$p_val_adj > 0])
  bulk_de_min = min(bulk_de$p_val_adj[bulk_de$p_val_adj > 0])
  sc_de %<>% 
    mutate(p_val_adj = ifelse(p_val_adj <= sc_de_min, sc_de_min, p_val_adj))
  bulk_de %<>% 
    mutate(p_val_adj = ifelse(p_val_adj <= bulk_de_min, bulk_de_min, p_val_adj))
  ## repeat for raw p-values
  sc_de_min = min(sc_de$p_val[sc_de$p_val > 0])
  bulk_de_min = min(bulk_de$p_val[bulk_de$p_val > 0])
  sc_de %<>% 
    mutate(p_val = ifelse(p_val <= sc_de_min, sc_de_min, p_val))
  bulk_de %<>%
    mutate(p_val = ifelse(p_val <= bulk_de_min, bulk_de_min, p_val))
  
  # filter to genes detected in both single-cell and bulk data
  genes = intersect(bulk_de$gene, sc_de$gene)
  sc_de %<>% filter(gene %in% genes) %>% arrange(gene)
  bulk_de %<>% filter(gene %in% genes) %>% arrange(gene)
  
  if (method == 'fcc') {
    # fold-change correlation
    genes = intersect(bulk_de$gene, sc_de$gene)
    bulk_de %<>% filter(gene %in% genes)
    sc_de %<>% filter(gene %in% genes)
    cor = cor(
      bulk_de %>%
        mutate(stat = sign(avg_logFC) * abs(test_statistic)) %>% 
        arrange(gene) %>%
        pull(stat),
      sc_de %>%
        mutate(stat = sign(avg_logFC) * abs(test_statistic)) %>% 
        arrange(gene) %>%
        pull(stat),
      method = cor_method
    )
  } else if (method == 'aucc') {
    # area under the concordance curve
    k = as.integer(k)
    ## rank in descending order first by p_val
    ## break ties by the abs() of the test_statistic
    vec1 = bulk_de %>%
      arrange(p_val, desc(abs(test_statistic))) %>%
      pull(gene) %>%
      head(k)
    vec2 = sc_de %>%
      arrange(p_val, desc(abs(test_statistic))) %>%
      pull(gene) %>%
      head(k)
    
    concordance_curve = map_dbl(seq_len(k), ~ {
      v1 = vec1[seq_len(.)]
      v2 = vec2[seq_len(.)]
      length(intersect(v1, v2))
    })
    denom = k * (k + 1) / 2
    aucc = sum(concordance_curve) / denom
    return(aucc)
  }
}
