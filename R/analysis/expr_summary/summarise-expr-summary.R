setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(data.table)
args = list(); source("R/functions/detect_system.R")

# list input files
summary_dir = file.path(base_dir, "analysis", "expr_summary")
summary_files = list.files(summary_dir, full.names = TRUE, pattern = "*gz")

# read them all
dats = map(summary_files, fread)

# combine
dat = do.call(rbind, dats)

# pick one comparison per dataset
confounds = readRDS("data/analysis/confounds/confounds.rds") %>%
  # fix a couple datasets
  mutate(comparison = ifelse(grepl("Schafflick|Der", dataset),
                             gsub("^.*_", "", dataset),
                             ifelse(grepl("Hagai", dataset), 
                                    paste0(gsub("^.*_", "", dataset), "|", 
                                           comparison),
                                    comparison)),
         dataset = ifelse(grepl("Schafflick|Der|Hagai", dataset),
                          gsub("_.*$", "", dataset), dataset))
# pick the comparison with the most cells
most_cells = confounds %>%
  filter(outcome == '# of cells') %>%
  group_by(dataset, comparison) %>%
  summarise(total_cells = sum(value),
            n_cell_types = n_distinct(cell_type)) %>%
  ungroup() %>% 
  group_by(dataset) %>%
  arrange(desc(total_cells), desc(n_cell_types)) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(dataset, comparison)
n_distinct(most_cells$dataset)

# filter to these comparisons
dat0 = dat %>%
  mutate(comparison = ifelse(grepl("Schafflick|Der", dataset),
                             gsub("^.*_", "", dataset),
                             ifelse(grepl("Hagai", dataset), 
                                    paste0(gsub("^.*_", "", dataset), "|", 
                                           comparison),
                                    comparison)),
         dataset = ifelse(grepl("Schafflick|Der|Hagai", dataset),
                          gsub("_.*$", "", dataset), dataset)) %>%
  inner_join(most_cells, by = c('dataset', 'comparison'))

# write
saveRDS(dat0, "data/analysis/expr_summary/expr_summary.rds")
