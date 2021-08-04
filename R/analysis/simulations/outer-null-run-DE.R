setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-nullrun-DE.R')
parser$add_argument('--allocation', type = 'character')
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions/datasets.R")
source("R/functions/submit_job.R")
source("R/functions/detect_system.R")

# list input files
input_dir = file.path(base_dir, "analysis", "simulations", "null", "objects")

grid = tidyr::crossing(
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
    n_reps == 6 & n_cells %in% c(100, 200, 500, 1000, 2000))

# write the raw array
grid_file = "sh/analysis/simulations/grids/null_run_DE.raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "analysis/simulations/null/DE")

# check which parameters are already complete
overwrite = F
grid0 = grid
if (!overwrite) {
  grid0 = grid %>%
    mutate(output_filename = paste0("GSE96583",
                                   "-n_cells=", n_cells,
                                   "-de_prob=", de_prob,
                                   "-de_loc=", de_loc,
                                   "-n_reps=", n_reps,
                                   "-sample_idx=", sample_idx,
                                    '-de_test=', de_test,
                                    '-shuffle_replicates=', shuffle_replicates,
                                    '.rds'),
          input_filename = paste0("GSE96583",
                                         "-n_cells=", n_cells,
                                         "-de_prob=", de_prob,
                                         "-de_loc=", de_loc,
                                         "-n_reps=", n_reps,
                                         "-sample_idx=", sample_idx,
                                          '.rds'),
           output_file = file.path(output_dir, output_filename),
           input_file = file.path(input_dir, input_filename),
           exists = file.exists(output_file)) %>%
    filter(!exists) %>%
    dplyr::select(-output_file, -output_filename, -exists) %>%
    # clean up grid
    dplyr::select(input_file, de_test, shuffle_replicates)
}

# run 5000 at a time
grid0 %<>% dplyr::slice(1:5000)

# write the grid that still needs to be run
write.table(grid0, "sh/analysis/simulations/grids/null_run_DE.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
sh_dir = "~/git/DE-analysis/sh/analysis/simulations"
script = file.path(sh_dir, "null_run_DE.sh")
submit_job(grid0, script, args$allocation, system)
