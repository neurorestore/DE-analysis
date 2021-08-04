# Run bulk DE analysis on a dataset.
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-run-DE.R')
parser$add_argument('--allocation', type = 'character')
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions/datasets.R")
source("R/functions/submit_job.R")
source("R/functions/detect_system.R")

# list input files
input_files = file.path(base_dir, paste0(bulk_datasets, '.rds'))
inputs = data.frame(input_file = input_files) %>%
  mutate(type = dirname(bulk_datasets))

# define tests to use
de_tests = c("bulk_DESeq2,test?LRT",
             "bulk_DESeq2,test?Wald",
             "bulk_limma,mode?voom",
             "bulk_limma,mode?trend",
             "bulk_edgeR,test?LRT",
             "bulk_edgeR,test?QLF")

# rep analysis grid over input files
grid = inputs %>%
  dplyr::slice(rep(1:n(), each = length(de_tests))) %>%
  mutate(de_test = rep(de_tests, nrow(inputs))) %>%
  # only do limma for proteomics and microarray
  filter(type != 'proteomics' | !grepl("DESeq2|edgeR", de_test))

# write the raw array
grid_file = "sh/analysis/run_DE/grids/run_bulk_DE.raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "analysis/run_bulk_DE")

# check which parameters are already complete
overwrite = F
grid0 = grid
if (!overwrite) {
  grid0 = grid %>%
    mutate(output_filename = paste0(basename(input_file) %>%
                                      gsub("\\.rds$", "", .),
                                    '-de_test=', de_test,
                                    '.rds'),
           output_file = file.path(output_dir, output_filename),
           exists = file.exists(output_file)) %>%
    filter(!exists) %>%
    dplyr::select(-output_file, -output_filename, -exists)
}

# write the grid that still needs to be run
write.table(grid0, "sh/analysis/run_DE/grids/run_bulk_DE.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
sh_dir = "~/git/DE-analysis/sh/analysis/run_DE"
script = file.path(sh_dir, "run_bulk_DE.sh")
submit_job(grid0, script, args$allocation, system)
