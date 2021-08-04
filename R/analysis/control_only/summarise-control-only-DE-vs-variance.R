setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(data.table)
args = list(); source("R/functions/detect_system.R")

# list input files
input_dir = file.path(base_dir, "analysis", "expr_summary", "control_only")
input_files = list.files(input_dir, pattern = '*\\.csv\\.gz', full.names = TRUE)

# we don't need to summarize all data, so let's filter here
meta = data.frame(filename = basename(input_files)) %>%
  mutate(idx = row_number()) %>%
  separate(filename, into = c('dataset', 'label'), sep = '-') %>%
  mutate_all(~ gsub("^.*=|\\.csv\\.gz", "", .)) %>%
  type_convert()
# filter to control groups in simple experiments only
keep = c('Goldfarbmuren2020' = 'never', ## lung from never smokers
         'Grubman2019' = 'Control', ## ALZ control brains
         'Hrvatin2018' = '0h', ## light-deprived mice
         'Huang2020' = 'control', ## colonic mucosa in healthy children
         'Kang2018' = 'ctrl', ## unstimulated PBMCs 
         'Mathys2019' = 'Control', ## ALZ control brains
         'Nagy2020' = 'Control', ## MDD control brains
         'Reyfman2020' = 'Control', ## healthy lungs
         'Rossi2019' = 'control', ## mice on a control diet
         'Sathyamurthy2018' = 'control', ## healthy mouse spinal cord
         'Smillie2019' = 'Healthy', ## healthy colon
         'Tran2019' = 'Ctrl', ## uninjured RGCs
         'Wilk2020' = 'Healthy', ## control PBMCs
         'Wu2017' = 'control' ## control mice
) %>%
  data.frame(dataset = names(.), label = .)

# filter metadata/files accordingly
meta0 = inner_join(meta, keep, by = c('dataset', 'label'))
input_files %<>% extract(meta0$idx)

# read all data
dats = map(input_files, fread)
# combine into a single data frame
dat = bind_rows(dats)

# last, we also need to load the DE results
DE = readRDS(file.path(base_dir, "analysis", "summary_data", 
                       "control_only.rds"))
n_DE = readRDS("data/analysis/control_only/n_DE_genes.rds")

## outcome 1: write mean delta-variance for each cell type in each dataset
delta_vars = dat %>%
  drop_na(pseudobulk_variance, shuffled_variance) %>%
  mutate(delta_variance = shuffled_variance - pseudobulk_variance) %>%
  group_by(dataset, label, cell_type) %>%
  summarise(mean_delta = mean(delta_variance)) %>%
  ungroup()
saveRDS(delta_vars, "data/analysis/control_only/delta_variance.rds")

## outcome 2: number of DE genes in each bin
bins = 10
xy0 = xy %>% 
  mutate(abs_delta_variance = abs(delta_variance))
bin_results = xy0 %>%
  # bin expression levels
  group_by(dataset, label, cell_type, de_test, shuffle_replicates) %>%
  arrange(abs_delta_variance) %>%
  mutate(bin = cut(row_number() / n(),
                   breaks = seq(0, bins) / bins),
         bin = as.integer(bin)) %>%
  ungroup() %>%
  # count DE genes in each bin
  group_by(dataset, label, cell_type, de_test, shuffle_replicates, bin) %>%
  summarise(genes = sum(p_val_adj < 0.05)) %>%
  ungroup()
saveRDS(bin_results, "data/analysis/control_only/genes_per_bin.rds")
