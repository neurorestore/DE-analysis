setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/functions/recode_colnames.R")
args = list(); source("R/functions/detect_system.R")

# list input files
input_dir = file.path(base_dir, "analysis", "control_only", "spatial")
input_files = list.files(input_dir, pattern = 'rds', full.names = TRUE)

meta = data.frame(filename = basename(input_files)) %>%
  mutate(idx = row_number()) %>%
  separate(filename, into = c('dataset', 'de_test', 'sample_idx',
                              'shuffle_replicates', 'label'), sep = '-') %>%
  mutate_all(~ gsub("^.*=|\\.rds", "", .)) %>%
  type_convert() %>%
  # remove superfluous columns
  dplyr::select(-sample_idx)

# read all data
dats = map(input_files, ~ readRDS(.x) %>%
             bind_rows() %>%
             # fix column names
             recode_colnames() %>%
             # fix p-values
             group_by(comparison, cell_type) %>%
             mutate(p_val_adj = p.adjust(p_val, method = 'BH')) %>%
             ungroup())

# combine into a single data frame
dat = meta %>%
  split(.$idx) %>%
  map2(dats, ~ cbind(.x, .y)) %>%
  bind_rows()
# reorder columns
dat0 = dat %>%
  dplyr::select(dataset, comparison, label, cell_type,
                de_test, shuffle_replicates,
                gene, p_val, test_statistic, p_val_adj) %>%
  # remove missing genes
  drop_na(p_val)

# save the full set of results
saveRDS(dat0, file.path(base_dir, "analysis", "summary_data",
                        "control_only_spatial.rds"))

# count the number of genes 
n_genes = dat0 %>%
  group_by(dataset, comparison, label, cell_type, de_test, 
           shuffle_replicates) %>%
  summarise(n = sum(p_val_adj < 0.05)) %>%
  ungroup()

# write # of DE genes
output_file = "data/analysis/control_only/spatial/n_DE_genes.rds"
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir)
saveRDS(n_genes, output_file)
