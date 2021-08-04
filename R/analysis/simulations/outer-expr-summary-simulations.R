setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-expr-summary-simulations.R')
parser$add_argument('--allocation', type = 'character')
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions/datasets.R")
source("R/functions/submit_job.R")
source("R/functions/detect_system.R")

# set up grid
# limit this experiment to n_reps=3
input_dir = file.path(base_dir, "analysis", "simulations", "null", "objects")
grid = tidyr::crossing(
  n_cells = c(100, 200, 500, 1000, 2000),
  de_prob = 0.5,
  de_loc = seq(0, 1, 0.1),
  n_reps = c(3, 4, 5, 10, 20) * 2,
  sample_idx = seq_len(10)
) %>%
  # vary only one of n_cells/n_reps
  filter(n_cells == 500 | n_reps == 6) %>%
  mutate(input_filename = paste0("GSE96583",
                                 "-n_cells=", n_cells,
                                 "-de_prob=", de_prob,
                                 "-de_loc=", de_loc,
                                 "-n_reps=", n_reps,
                                 "-sample_idx=", sample_idx,
                                 '.rds'),
         input_file = file.path(input_dir, input_filename))

# define output directory where results are stored
output_dir = file.path(base_dir, "analysis/simulations/null/expr_summary")

# check which parameters are already complete
overwrite = F
grid0 = grid
if (!overwrite) {
  grid0 = grid %>%
    mutate(output_filename = basename(input_file),
           output_file = file.path(output_dir, output_filename),
           exists = file.exists(output_file)) %>%
    filter(!exists) %>%
    dplyr::select(input_file)
}

# write the grid that still needs to be run
grid_file = "sh/analysis/simulations/grids/expr_summary.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid0, grid_file, quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
sh_dir = "~/git/DE-analysis/sh/analysis/simulations"
script = file.path(sh_dir, "expr_summary_simulations.sh")
submit_job(grid0, script, args$allocation, system)
