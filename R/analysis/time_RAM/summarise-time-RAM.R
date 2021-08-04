setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
args = list(); source("R/functions/detect_system.R")

# set up input
input_dir = file.path(base_dir, "analysis/run_DE")
input_files = list.files(input_dir, full.names = T, pattern ='*\\.rds$')

# extract walltime and RAM from all input files
dats = map(input_files, ~ {
  dat = readRDS(.)
  map(dat, ~ distinct(., runtime, mem_usage)) %>%
    bind_rows(.id = 'comparison')
}) %>%
  setNames(basename(input_files))

# combine all results
res = dats %>%
  bind_rows(.id = 'filename') %>%
  separate(filename, into = c('dataset', 'de_test', 'shuffle_replicates'),
           sep = '-') %>%
  mutate_at(vars(de_test, shuffle_replicates), ~ gsub("^.*=|\\.rds$", "", .)) %>%
  type_convert()

# print summary
res %>%
  group_by(de_test) %>%
  summarise(mean_ram = mean(mem_usage / 1e3),
            mean_time = mean(runtime / 60)) %>%
  arrange(desc(mean_time))

# save results
output_file = "data/analysis/run_DE/time_RAM.rds"
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)
saveRDS(res, output_file)
