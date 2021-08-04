# Tally the number of DE genes within bins of genes grouped by delta-variance.
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/functions/recode_colnames.R")
args = list(); source("R/functions/detect_system.R")

# set up grid
# limit this experiment to n_reps=3, n_cells=500
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
  n_cells = 500,
  de_prob = 0.5,
  de_loc = seq(0, 1, 0.1),
  n_reps = 6,
  sample_idx = seq_len(10),
  shuffle_replicates = c("NO", "YES")
) %>%
  filter(grepl("pseudo|mixed", de_test) | shuffle_replicates == 'NO') %>%
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

# also set up the expr_summary files
summary_dir = file.path(base_dir, "analysis", "simulations", "null",
                        "expr_summary")
summary_files = tidyr::crossing(
  n_cells = 500,
  de_prob = 0.5,
  de_loc = seq(0, 1, 0.1),
  n_reps = 6,
  sample_idx = seq_len(10)
) %>% mutate(summary_filename = paste0("GSE96583",
                                       "-n_cells=", n_cells,
                                       "-de_prob=", de_prob,
                                       "-de_loc=", de_loc,
                                       "-n_reps=", n_reps,
                                       "-sample_idx=", sample_idx,
                                       '.rds'),
             summary_file = file.path(summary_dir, summary_filename)) %>%
  pull(summary_file)

# read all data
dats = map(input_files, ~ readRDS(.x) %>%
             # fix column names
             recode_colnames() %>%
             # fix p-values
             group_by(cell_type) %>%
             mutate(p_val_adj = p.adjust(p_val, method = 'BH')) %>%
             ungroup()
) %>% setNames(basename(input_files))

# read all summary files
summary_dats = map(summary_files, readRDS) %>%
  setNames(basename(summary_files))

# combine all the data and join the two sources of data together
DE = bind_rows(dats, .id = 'filename') %>%
  separate(filename, into = c('dataset', 'n_cells', 'de_prob', 'de_loc',
                              'n_reps', 'sample_idx', 'de_test', 
                              'shuffle_replicates'), sep = '-') %>%
  mutate_all(~ gsub("^.*=|\\.rds$", "", .)) %>%
  type_convert() %>%
  # remove some columns
  dplyr::select(-test, -runtime, -mem_usage, -baseMean, -lfcSE, -logCPM,
                -AveExpr, -B, -used_voom)
summary = bind_rows(summary_dats, .id = 'filename') %>%
  separate(filename, into = c('dataset', 'n_cells', 'de_prob', 'de_loc',
                              'n_reps', 'sample_idx'), sep = '-') %>%
  mutate_all(~ gsub("^.*=|\\.rds$", "", .)) %>%
  type_convert()
dat = left_join(DE, summary, by = c('dataset', 'n_cells', 'de_prob', 'de_loc',
                                    'n_reps', 'sample_idx', 'gene'))

# save the complete dataset
saveRDS(dat, file.path(base_dir, "analysis/simulations/null/expr_summary.rds"))

# now, calculate number of DE genes per bin
bins = 10
dat0 = dat %>% 
  mutate(delta_variance = shuffled_variance - pseudobulk_variance,
         abs_delta_variance = abs(delta_variance))
bin_results = dat0 %>%
  # bin expression levels
  group_by(dataset, n_cells, de_prob, de_loc, n_reps, sample_idx, de_test,
           shuffle_replicates) %>%
  arrange(abs_delta_variance) %>%
  mutate(bin = cut(row_number() / n(),
                   breaks = seq(0, bins) / bins),
         bin = as.integer(bin)) %>%
  ungroup() %>%
  # count DE genes in each bin
  group_by(dataset, n_cells, de_prob, de_loc, n_reps, sample_idx, de_test,
           shuffle_replicates, bin) %>%
  summarise(genes = sum(p_val_adj < 0.05)) %>%
  ungroup() 

# save results
output_file = "data/analysis/simulations/null-genes-per-bin.rds"
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)
saveRDS(bin_results, output_file)
