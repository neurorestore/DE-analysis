setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/functions/datasets.R")
args = list(); source("R/functions/detect_system.R")

# list input files
input_dir = file.path(base_dir, "analysis", "delta_variance")
input_files = file.path(input_dir, paste0(datasets, '.rds'))

# read all files
dats = map(input_files, readRDS)
dat = bind_rows(dats) %>% type_convert()

# save results
output_file = "data/analysis/delta_variance/delta_variance.rds"
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
saveRDS(dat, output_file)
