# Count the total number of DE genes in each control-only experiment.
setwd("~/git/DE_analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
source("R/functions/recode_colnames.R")
args = list(); source("R/functions/detect_system.R")

# list input files
input_dir = file.path(base_dir, "analysis", "control_only")
input_files = list.files(input_dir, pattern = 'rds', full.names = TRUE)

# we don't need to summarize all data, so let's filter here
meta = data.frame(filename = basename(input_files)) %>%
  mutate(idx = row_number()) %>%
  separate(filename, into = c('dataset', 'de_test', 'sample_idx', 
                              'shuffle_replicates', 'label', 'comparison'), 
           sep = '-') %>%
  mutate_all(~ gsub("^.*=|\\.rds", "", .)) %>%
  type_convert() %>%
  # remove superfluous columns
  dplyr::select(-sample_idx)
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
dats = map(input_files, ~ readRDS(.x) %>%
             # fix column names
             recode_colnames() %>%
             # fix p-values
             group_by(cell_type) %>%
             mutate(p_val_adj = p.adjust(p_val, method = 'BH')) %>%
             ungroup()
)

# combine into a single data frame
dat = meta0 %>%
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
# remove duplicated data from multiple comparisons
# these are irrelevant since DE takes place within controls
dat0 %<>%
  group_by(dataset) %>%
  filter(comparison == first(comparison)) %>%
  ungroup()

# save the full set of results
saveRDS(dat0, file.path(base_dir, "analysis", "summary_data",
                        "control_only.rds"))

# count the number of genes 
n_genes = dat0 %>%
  group_by(dataset, comparison, label, cell_type, de_test, 
           shuffle_replicates) %>%
  summarise(n = sum(p_val_adj < 0.05)) %>%
  ungroup()

# write # of DE genes
output_file = "data/analysis/control_only/n_DE_genes.rds"
output_dir = dirname(output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir)
saveRDS(n_genes, output_file)
