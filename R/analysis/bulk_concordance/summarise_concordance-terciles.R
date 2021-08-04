setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
args = list(); source("R/functions/detect_system.R")

# set up input
input_dir = file.path(base_dir, "analysis/run_concordance_terciles")
input_files = list.files(input_dir, full.names = T, pattern ='*\\.rds$')

# read all input files
dats = map(input_files, readRDS) %>%
  setNames(basename(input_files))

# combine into a single file
dat = dats %>%
  bind_rows(.id = 'comparison') %>%
  separate(comparison, c("sc", "bulk"), "\\|") %>%
  separate(sc, c("sc_dataset", "sc_test", "shuffle_replicates"), "-") %>%
  separate(bulk, c("bulk_dataset", "bulk_test", "n_bins"), "-") %>%
  mutate_at(vars(sc_test, bulk_test, shuffle_replicates, n_bins), function(x) 
    gsub(".*=|.rds", "", x))

# save results
output_file = "data/analysis/bulk_concordance/concordance_terciles.rds"
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)
saveRDS(dat, output_file)
