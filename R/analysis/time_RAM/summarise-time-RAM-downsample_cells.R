# Summarize the time and RAM usage of the default DE analyses in datasets
# downsampled to a fixed number of cells (used to test mixed models).
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
args = list(); source('R/functions/detect_system.R')

# set up input directory
input_dir = file.path(base_dir, "analysis/downsample_cells/DE")
input_files = list.files(input_dir, full.names = T, pattern ='*\\.rds$') %>%
  extract(grepl("Angelidis|Hagai|CanoGamez|Reyfman", .))

# extract walltime and RAM from all input files
dats = map(input_files, ~ {
  print(.)
  dat = readRDS(.) %>%
    # remove empty data frames
    extract(map_int(., nrow) > 0)
  map(dat, ~ distinct(., runtime, mem_usage)) %>%
    bind_rows(.id = 'comparison')
}) %>%
  setNames(basename(input_files))

# combine all results
res = dats %>%
  bind_rows(.id = 'filename') %>%
  separate(filename, into = c('dataset', 'de_test', 'n_cells', 'sample_idx'),
           sep = '-') %>%
  mutate_at(vars(de_test, n_cells, sample_idx), 
            ~ gsub("^.*=|\\.rds$", "", .)) %>%
  type_convert()

# print summary
res %>%
  group_by(de_test) %>%
  summarise(mean_ram = mean(mem_usage / 1e3),
            mean_time = mean(runtime / 60)) %>%
  arrange(desc(mean_time))

# save results
output_file = "data/analysis/downsample_cells/time_RAM.rds"
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)
saveRDS(res, output_file)
