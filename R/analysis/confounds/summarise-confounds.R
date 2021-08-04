setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
args = list(); source("R/functions/detect_system.R")

# list files
input_dir = file.path(base_dir, "analysis/confounds")
input_files = list.files(input_dir, full.names = T, pattern = '*\\.txt\\.gz$')

# read these all
dats = map(input_files, read.csv) 
dat = do.call(rbind, dats)

# save results
output_file = "data/analysis/confounds/confounds.rds"
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)
saveRDS(dat, output_file)
