# Run gene set enrichment analysis (GSEA) on single-cell DE results.
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-run-GSEA.R')
parser$add_argument('--allocation', type = 'character')
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions/submit_job.R")
source("R/functions/detect_system.R")

# manually set up the input single-cell datasets
sc_datasets = c(paste0('Hagai2018_', c('rat', 'rabbit', 'mouse', 'pig')),
                'CanoGamez2020',
                'Angelidis2019',
                'Reyfman2020')

# establish analysis grid
opts = list(
  dataset = sc_datasets,
  de_test = c(
    "wilcox", "bimod", "t", "negbinom", "poisson", "LR", "MAST",
    "pseudobulk_DESeq2,test?LRT", "pseudobulk_DESeq2,test?Wald",
    "pseudobulk_limma,mode?voom", "pseudobulk_limma,mode?trend",
    "pseudobulk_edgeR,test?QLF", "pseudobulk_edgeR,test?LRT",
    "mixed_lm")
)
sc_grid = do.call(tidyr::crossing, opts) %>%
  mutate(input_dir = file.path(base_dir, "analysis", "run_DE"),
         output_dir = file.path(base_dir, "analysis", "run_GSEA",
                                "single_cell"))

# now, do the same for bulk grid
bulk_datasets = c(paste0('Hagai2018_', c('mouse', 'pig', 'rat', 'rabbit')),
                  'Angelidis2019_facsepi',
                  'Angelidis2019_facsmac',
                  'CanoGamez2020',
                  'Reyfman2020_alvmac',
                  'Reyfman2020_AT2')
opts = list(
  dataset = bulk_datasets,
  de_test = c("bulk_DESeq2,test?LRT",
              "bulk_DESeq2,test?Wald",
              "bulk_limma,mode?voom",
              "bulk_limma,mode?trend",
              "bulk_edgeR,test?LRT",
              "bulk_edgeR,test?QLF")
)
bulk_grid = do.call(tidyr::crossing, opts) %>%
  mutate(input_dir = file.path(base_dir, "analysis", "run_bulk_DE"),
         output_dir = file.path(base_dir, "analysis", "run_GSEA", "bulk"))

# combine grids
grid = bind_rows(sc_grid, bulk_grid)

# write the raw array
grid_file = "sh/analysis/run_GSEA/grids/run_GSEA.raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# check which parameters are already complete
overwrite = F
grid0 = grid
if (!overwrite) {
  grid0 = grid %>%
    # for Reyfman2020, recode bulk test
    mutate(de_test = ifelse(grepl("Reyfman2020_", dataset),
                            'bulk_DESeq2', de_test)) %>%
    distinct() %>%
    mutate(input_file = file.path(input_dir, paste0(dataset, 
                                                    '-de_test=', de_test, 
                                                    '.rds')),
           output_file = file.path(output_dir,  paste0(dataset,
                                                       '-de_test=', de_test,
                                                       '.rds')),
           exists = file.exists(output_file)) %>%
    filter(!exists) %>%
    dplyr::select(dataset, de_test, input_file, output_file)
}

# subset grid, if needed
if (nrow(grid0) >= 10000) {
  grid0 %<>% dplyr::slice(1:9900) ## allow for some other running jobs or sh
}

# write the grid that still needs to be run
write.table(grid0, "sh/analysis/run_GSEA/grids/run_GSEA.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
sh_dir = "~/git/DE-analysis/sh/analysis/run_GSEA"
script = file.path(sh_dir, "run_GSEA.sh")
submit_job(grid0, script, args$allocation, system)
