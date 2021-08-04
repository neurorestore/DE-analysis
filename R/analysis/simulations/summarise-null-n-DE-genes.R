setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)

# set up input directory
source("R/functions/detect_system.R")

# first, summarise the effect of n_reps, at n_cells == 500
input_dir = file.path(base_dir, "analysis", "simulations", "null", "DE")
input_files = tidyr::crossing(
  de_test = c(
    ## single-cell methods, implemented in Seurat
    "wilcox", "bimod", "t", "negbinom", "poisson", "LR", "MAST",
    # pseudobulk methods
    "pseudobulk_DESeq2,test?LRT", "pseudobulk_DESeq2,test?Wald",
    "pseudobulk_limma,mode?voom", "pseudobulk_limma,mode?trend",
    "pseudobulk_edgeR,test?QLF",
    "pseudobulk_edgeR,test?LRT",
    # mixed model, implemented in Seurat
    "mixed_lm"
  ),
  n_cells = c(100, 200, 500, 1000, 2000),
  de_prob = 0.5,
  de_loc = seq(0, 1, 0.1),
  n_reps = 2 * c(3, 4, 5, 10, 20),
  sample_idx = seq_len(10),
  shuffle_replicates = c("NO", "YES")
) %>%
  filter(grepl("pseudo|mixed", de_test) | shuffle_replicates == 'NO') %>%
  filter(n_cells == 500 & n_reps %in% c(2 * c(3, 4, 5, 10, 20)) |
           n_reps == 6 & n_cells %in% c(100, 200, 500, 1000, 2000)) %>%
  mutate(input_filename = paste0("GSE96583",
                                 "-n_cells=", n_cells,
                                 "-de_prob=", de_prob,
                                 "-de_loc=", de_loc,
                                 "-n_reps=", n_reps,
                                 "-sample_idx=", sample_idx,
                                 '-de_test=', de_test,
                                 '-shuffle_replicates=', shuffle_replicates,
                                 '.rds'),
         input_file = file.path(input_dir, input_filename)) %>%
  pull(input_file)

# read all input files
dats = map(input_files, readRDS) %>%
  setNames(basename(input_files))

# calculate # of DE genes
n_DE = dats %>%
  map(~ {
    DE = bind_rows(., .id = 'comparison')
    # fix column names
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
    # re-calculate adjusted p-values using BH correction (Seurat does Bonferroni)
    DE %<>%
      group_by(comparison, cell_type) %>%
      mutate(p_val_adj = p.adjust(p_val, method = 'BH')) %>%
      ungroup()
    # combine results
    DE %>%
      group_by(comparison, test, cell_type) %>%
      summarise(n_1 = sum(p_val_adj < 0.01, na.rm = T),
                n_5 = sum(p_val_adj < 0.05, na.rm = T),
                n_10 = sum(p_val_adj < 0.1, na.rm = T)) %>%
      ungroup() %>%
      gather('fdr', 'n_genes', n_1:n_10) %>%
      mutate(fdr = paste0(gsub("^.*_", "", fdr), "%"))
  }) %>%
  bind_rows(.id = 'filename') %>%
  separate(filename, into = c(
    'dataset',
    'n_cells',
    'de_prob',
    'de_loc',
    'n_reps',
    'sample_idx',
    'de_test',
    'shuffle_replicates'),
    sep = '-') %>%
  mutate_at(vars(n_cells, de_prob, de_loc, n_reps, sample_idx,
                 de_test, shuffle_replicates), ~ gsub("^.*=|\\.rds$", "", .))

# save results
output_file = "data/analysis/simulations/null-n-DE-genes.rds"
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)
saveRDS(n_DE, output_file)
