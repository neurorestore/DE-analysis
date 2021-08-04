setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
args = list(); source("R/functions/detect_system.R")

input_dir = file.path(base_dir, "analysis/extract_FPs")
input_files = list.files(input_dir, full.names = T, pattern ='*\\.rds$')

# read all input files
dats = map(input_files, readRDS) %>%
  setNames(basename(input_files))

# extract FPs/FNs separately
FPs = map(dats, 'FPs')
FNs = map(dats, 'FNs')

# combine into a single object
FP = FPs %>%
  map(~ map_dfc(., as.character)) %>%
  bind_rows(.id = 'filename') %>%
  type_convert() %>%
  # extract missing info from filename
  separate(filename, into = c('dataset', 'sc_test', 'shuffle_replicate', 
                              'bulk_test'), sep = '-') %>%
  mutate_at(vars(sc_test, shuffle_replicate, bulk_test), function(x) 
    gsub(".*=|.rds", "", x))
FN = FNs %>%
  map(~ map_dfc(., as.character)) %>%
  bind_rows(.id = 'filename') %>%
  type_convert() %>%
  # extract missing info from filename
  separate(filename, into = c('dataset', 'sc_test', 'shuffle_replicate',
                              'bulk_test'), sep = '-') %>%
  mutate_at(vars(sc_test, shuffle_replicate, bulk_test), function(x) 
    gsub(".*=|.rds", "", x))

# create output
dat = list(FPs = FP, FNs = FN)

# save results
output_file = "data/analysis/extract_FPs/extract_FPs.rds"
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)
saveRDS(dat, output_file)
