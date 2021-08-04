setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
args = list(); source("R/functions/detect_system.R")

# set up input directory
input_dir = file.path(base_dir, "analysis/run_spike_in_DE")
input_files = list.files(input_dir, full.names = T, pattern ='*\\.rds$')

# read all input files
dats = map(input_files, readRDS) %>%
  setNames(basename(input_files))

## function to fix column names
clean_cols = function(df) {
  # fix column names
  colnames(df) %<>%
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
  return(df)  
}

# summarise, keeping only ERCCs
sum = dats %>%
  map(extract2, 1) %>%
  map(clean_cols) %>%
  bind_rows(.id = 'dataset') %>%
  filter(grepl("^ERCC-", gene)) %>%
  dplyr::select(dataset, cell_type, gene, p_val, p_val_adj, test,
                test_statistic, avg_logFC)

# combine this with gene level summary statistics
expr_summary = read.csv(file.path(base_dir, "analysis/expr_summary",
                                  "Hagai2018_plate.txt.gz")) %>%
  filter(gene %in% sum$gene)
sum %<>% left_join(expr_summary, by = 'gene') %>%
  dplyr::select(-dataset.y, -cell_type.x) %>%
  dplyr::rename(cell_type = cell_type.y) %>%
  separate(dataset.x, into = c('dataset', 'de_test', 'shuffle_replicates'), 
           sep = '-') %>%
  mutate_at(vars(de_test, shuffle_replicates), ~ gsub("^.*=|\\.rds$", "", .))

# save results
output_file = "data/analysis/run_spike_in_DE/spike_in_summary.rds"
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)
saveRDS(sum, output_file)
