setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-expr-summary.R')
parser$add_argument('--allocation', type = 'character')
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions/datasets.R")
source("R/functions/submit_job.R")
source("R/functions/detect_system.R")

# list input files
input_dir = file.path(base_dir, "rnaseq", "seurat")
input_files = file.path(input_dir, paste0(datasets, '.rds'))
# add Hagai plate data into this
input_files %<>% c(file.path(input_dir, "Hagai2018_plate.rds"))
# grid is simply the list of input files
grid = data.frame(input_file = input_files)

# define output directory where results are stored
output_dir = file.path(base_dir, "analysis/expr_summary")

# check which parameters are already complete
overwrite = F
grid0 = grid
if (!overwrite) {
  grid0 = grid %>%
    mutate(output_filename = paste0(basename(input_file) %>%
                                      gsub("\\.rds$", "", .), '.txt.gz'),
           output_file = file.path(output_dir, output_filename),
           exists = file.exists(output_file)) %>%
    filter(!exists) %>%
    dplyr::select(-output_file, -output_filename, -exists)
}

# write the grid that still needs to be run
grid_file = "sh/analysis/expr_summary/grids/expr_summary.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid0, grid_file, quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
sh_dir = "~/git/DE-analysis/sh/analysis/expr_summary"
script = file.path(sh_dir, "expr_summary.sh")
submit_job(grid0, script, args$allocation, system)
