setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
args = list(); source("R/functions/detect_system.R")

input_dir = file.path(base_dir, "analysis", "run_GSEA", "concordance")
input_files = list.files(input_dir, full.names = T, pattern ='*\\.rds$')

# read all input files
dats = map(input_files, readRDS) %>%
  setNames(basename(input_files))

# combine into a single file
dat = dats %>%
  map(~ map_dfc(., as.character)) %>%
  bind_rows(.id = 'comparison') %>%
  type_convert() %>%
  separate(comparison, into = c('bulk_dataset', 'sc_test', 'bulk_test'),
           sep = '-') %>%
  mutate_at(vars(sc_test, bulk_test), function(x) gsub(".*=|.rds", "", x))

# save results
output_file = "data/analysis/run_GSEA/GSEA_concordance.rds"
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)
saveRDS(dat, output_file)
