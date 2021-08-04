setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
args = list(); source("R/functions/detect_system.R")

# set up input directory
input_dir = file.path(base_dir, "analysis/downsample_cells/concordance")
input_files = list.files(input_dir, full.names = T, pattern = '*\\.rds$')

# read all input files
dats = map(input_files, readRDS) %>%
  setNames(basename(input_files))

# combine into a single file
dat = dats %>%
  bind_rows(.id = 'comparison') %>%
  separate(comparison, c("sc", "bulk"), "\\|") %>%
  separate(sc, c("dataset", "de_test", "n_cells", "sample_idx"), sep = "-") %>%
  mutate(de_test = gsub("!", "\\|", de_test)) %>%
  separate(bulk, c("bulk_dataset", "bulk_test"), "-") %>%
  mutate_at(vars(dataset, de_test, n_cells, sample_idx, bulk_test),
            ~ gsub(".*=|.rds", "", .)) %>%
  # remove old bulk test framework
  filter(!bulk_test %in% c("bulk_limma", "bulk_DESeq2", "bulk_edgeR") |
           # ... but keep published Reyfman analysis
           (bulk_test == 'bulk_DESeq2' & grepl("Reyfman", dataset)))

# save results
output_file = "data/analysis/downsample_cells/concordance_summary.rds"
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)
saveRDS(dat, output_file)
