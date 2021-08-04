setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-expr-summary-control-only.R')
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
inputs = data.frame(input_file = input_files)
grid = data.frame(input_file = input_files, dataset = datasets)

# load in number of replicates for each dataset
reps = readRDS("data/analysis/confounds/replicates.rds")

# grab the conditions where we have enough replicates
keep = reps %>%
  group_by(dataset, label, comparison) %>%
  summarise(n_reps = n()) %>%
  ungroup %>%
  filter(n_reps >= 6) %>%
  # ignore duplicate comparisons
  distinct(dataset, label, n_reps)

# add this into the grid
grid %<>%
  left_join(keep, by = 'dataset') %>%
  # drop the ones we aren't keeping
  drop_na() %>%
  dplyr::select(-n_reps)

# define output directory where results are stored
output_dir = file.path(base_dir, "analysis/expr_summary/control_only")

# check which parameters are already complete
overwrite = F
grid0 = grid
if (!overwrite) {
  grid0 = grid %>%
    mutate(output_filename = paste0(dataset, '-label=', label, '.csv.gz'),
           output_file = file.path(output_dir, output_filename),
           exists = file.exists(output_file)) %>%
    filter(!exists) %>%
    dplyr::select(-output_file, -output_filename, -exists,
                  -dataset)
}

# write the grid that still needs to be run
grid_file = "sh/analysis/control_only/grids/expr_summary.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid0, grid_file, quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
sh_dir = "~/git/DE-analysis/sh/analysis/control_only"
script = file.path(sh_dir, "expr_summary_control_only.sh")
submit_job(grid0, script, args$allocation, system)
