# Calculate the concordance GSEA results from matching single-cell and bulk DE.
setwd("~/git/DE-analysis")
options(stringsAsFactors = FALSE)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-GSEA-concordance.R')
parser$add_argument('--allocation', type = 'character')
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions/datasets.R")
source("R/functions/submit_job.R")
source("R/functions/detect_system.R")

# set up grid
opts = list(
  sc_dataset = c(paste0('Hagai2018_', c('mouse', 'pig', 'rat', 'rabbit')),
                 'Angelidis2019',
                 'CanoGamez2020', 
                 'Reyfman2020'),
  bulk_dataset = c(paste0('Hagai2018_', c('mouse', 'pig', 'rat', 'rabbit')),
                   'Angelidis2019_facsepi',
                   'Angelidis2019_facsmac',
                   'CanoGamez2020',
                   'Reyfman2020_alvmac',
                   'Reyfman2020_AT2'),
  sc_test = c("wilcox", "bimod", "t", "negbinom", "poisson", "LR", "MAST",
              "pseudobulk_DESeq2,test?LRT", "pseudobulk_DESeq2,test?Wald",
              "pseudobulk_limma,mode?voom", "pseudobulk_limma,mode?trend",
              "pseudobulk_edgeR,test?QLF", "pseudobulk_edgeR,test?LRT",
              "mixed_lm"),
  bulk_test = c("bulk_DESeq2,test?LRT",
                "bulk_DESeq2,test?Wald",
                "bulk_limma,mode?voom",
                "bulk_limma,mode?trend",
                "bulk_edgeR,test?LRT",
                "bulk_edgeR,test?QLF")
)
grid = do.call(tidyr::crossing, opts) %>%
  # matching datasets
  extract(map2_lgl(.$sc_dataset, .$bulk_dataset, ~ grepl(.x, .y)), )

# write the raw array
grid_file = "sh/analysis/run_GSEA/grids/GSEA_concordance.raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "analysis", "run_GSEA", "concordance")

# now, check for which parameters are already complete
overwrite = F
grid0 = grid
if (overwrite == F) {
  grid0 = grid %>%
    # for Reyfman2020, recode bulk test
    mutate(bulk_test = ifelse(grepl("Reyfman2020_", bulk_dataset),
                              'bulk_DESeq2', bulk_test)) %>%
    distinct() %>%
    # set up single-cell DE, bulk DE, and expr summary filepaths
    mutate(sc_dir = file.path(base_dir, 'analysis', 'run_GSEA', 'single_cell'),
           sc_filename = paste0(sc_dataset, '-de_test=', sc_test, '.rds'),
           sc_file = file.path(sc_dir, sc_filename),
           bulk_dir = file.path(base_dir, 'analysis', 'run_GSEA', 'bulk'),
           bulk_filename = paste0(bulk_dataset, '-de_test=', bulk_test, '.rds'),
           bulk_file = file.path(bulk_dir, bulk_filename)) %>%
    # set up output filepath
    mutate(output_filename = paste0(bulk_dataset, 
                                    '-sc_test=', sc_test,
                                    '-bulk_test=', bulk_test, 
                                    '.rds'),
           output_file = file.path(output_dir, output_filename),
           exists = file.exists(output_file),
           idx = row_number()) %>%
    # drop files that exist
    filter(!exists) %>%
    # keep only parameters and I/O
    dplyr::select(bulk_dataset, sc_file, bulk_file, output_file)
}

# write the grid that still needs to be run
write.table(grid0, "sh/analysis/run_GSEA/grids/GSEA_concordance.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
script = "~/git/DE-analysis/sh/analysis/run_GSEA/GSEA_concordance.sh"
submit_job(grid0, script, args$allocation, system)
