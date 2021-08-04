setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer_write_simulation-objects_null.R')
parser$add_argument('--allocation', type = 'character')
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions/datasets.R")
source("R/functions/submit_job.R")
source("R/functions/detect_system.R")

# define the output directory
output_dir = file.path(base_dir, "analysis/simulations/null/objects")

grid = tidyr::crossing(
  n_cells = c(100, 200, 500, 1000, 2000),
  de_prob = 0.5,
  de_loc = seq(0, 1, 0.1),
  n_reps = c(6, 8, 10, 20, 40),
  sample_idx = seq_len(10)
)

# write the raw array
grid_file = "sh/analysis/simulations/grids/write_null_objects.raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "analysis/simulations/null/objects")

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
                                    '.rds'),
           output_file = file.path(output_dir, output_filename),
           exists = file.exists(output_file)) %>%
    filter(!exists) %>%
    dplyr::select(-output_file, -output_filename, -exists)
}

# write the grid that still needs to be run
write.table(grid0, "sh/analysis/simulations/grids/write_null_objects.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
sh_dir = "~/git/DE-analysis/sh/analysis/simulations"
script = file.path(sh_dir, "write_null_objects.sh")
submit_job(grid0, script, args$allocation, system)
