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

# list bulk input files
bulk_files = list.files(file.path(base_dir, "analysis/run_bulk_DE"))
bulk_inputs = data.frame(bulk_file = bulk_files) %>%
  mutate(label = gsub("_.*|-.*", "", bulk_file)) %>%
  # manual fix for the Hagai datasets
  mutate(label = ifelse(label == 'Hagai2018', gsub("-.*", "", bulk_file),
                        label)) %>%
  # manual fix for the CanoGamez proteomics
  mutate(label = fct_recode(label,
                            "CanoGamez2020" = "CanoGamez2020:proteomics")) %>%
  # restore the entire filepath
  mutate(bulk_file = file.path(base_dir, 'analysis/run_bulk_DE', bulk_file))

# get single-cell comparison files
sc_files = list.files(file.path(base_dir, "analysis/run_DE"))
sc_inputs = data.frame(sc_file = sc_files) %>%
  mutate(label = gsub("_.*|-.*", "", sc_file)) %>%
  # manual fix for the Hagai datasets
  mutate(label = ifelse(label == 'Hagai2018', gsub("-.*", "", sc_file), label))

# rep analysis grid over input files
grid = bulk_inputs %>%
  left_join(sc_inputs) %>%
  drop_na()

# write the raw array
grid_file = "sh/analysis/run_DE/grids/run_concordance.raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "analysis/run_concordance")

# check which parameters are already complete
overwrite = F
grid0 = grid
if (!overwrite) {
  grid0 = grid %>%
    mutate(output_filename = paste0(basename(sc_file) %>%
                                      gsub("\\.rds$", "", .),
                                    "|",
                                    basename(bulk_file) %>%
                                      gsub("\\.rds$", "", .),
                                    '.rds'),
           output_file = file.path(output_dir, output_filename),
           exists = file.exists(output_file)) %>%
    filter(!exists) %>%
    dplyr::select(-output_file, -output_filename, -exists)
}

# limit to the 'gold standard' datasets
grid0 %<>% filter(grepl("Hagai|CanoGamez|Reyfman|Angelidis", label))

# write the grid that still needs to be run
write.table(grid0, "sh/analysis/run_DE/grids/run_concordance.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
sh_dir = "~/git/DE-analysis/sh/analysis/run_DE"
script = file.path(sh_dir, "run_concordance.sh")
submit_job(grid0, script, args$allocation, system)
