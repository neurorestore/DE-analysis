# Analyze the relationships between mean expression, expression variance, and
# delta-variance in all 46 datasets.
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(broom)
library(ppcor)

# read expr_summary data
dat = readRDS("data/analysis/expr_summary/expr_summary.rds") 

# for each dataset, calculate correlations between:
## mean and delta-variance
cors1 = dat %>%
  # first, correlate within cell types
  mutate(delta = shuffled_variance - pseudobulk_variance) %>%
  drop_na(mean, delta) %>%
  group_by(dataset, comparison, cell_type) %>%
  do(tidy(cor.test(.$mean, .$delta, method = 'p', use = 'p'))) %>%
  ungroup() %>%
  filter(is.finite(estimate)) %>%
  # next, average over cell types
  group_by(dataset, comparison) %>%
  summarise(mean_cor = mean(estimate, na.rm = TRUE)) %>%
  ungroup()
  
## variance and delta-variance
cors2 = dat %>%
  # first, correlate within cell types
  mutate(delta = shuffled_variance - pseudobulk_variance) %>%
  drop_na(pseudobulk_variance, delta) %>%
  group_by(dataset, comparison, cell_type) %>%
  do(tidy(cor.test(.$pseudobulk_variance, .$delta, method = 'p', use = 'p'))) %>%
  ungroup() %>%
  filter(is.finite(estimate)) %>%
  # next, average over cell types
  group_by(dataset, comparison) %>%
  summarise(mean_cor = mean(estimate, na.rm = TRUE)) %>%
  ungroup()

# now do partial correlations between:
## delta-variance and variance, controlling for mean
pcors1 = dat %>%
  # first, correlate within cell types
  mutate(delta = shuffled_variance - pseudobulk_variance) %>%
  drop_na(pseudobulk_variance, mean, delta) %>%
  group_by(dataset, comparison, cell_type) %>%
  mutate(partial_cor = pcor.test(pseudobulk_variance, delta, mean)$estimate) %>%
  ungroup() %>%
  filter(is.finite(partial_cor)) %>%
  # next, average over cell types
  group_by(dataset, comparison) %>%
  summarise(mean_cor = mean(partial_cor, na.rm = TRUE)) %>%
  ungroup()

## delta-variance and mean, controlling for variance
pcors2 = dat %>%
  # first, correlate within cell types
  mutate(delta = shuffled_variance - pseudobulk_variance) %>%
  drop_na(pseudobulk_variance, mean, delta) %>%
  group_by(dataset, comparison, cell_type) %>%
  mutate(partial_cor = pcor.test(mean, delta, pseudobulk_variance)$estimate) %>%
  ungroup() %>%
  filter(is.finite(partial_cor)) %>%
  # next, average over cell types
  group_by(dataset, comparison) %>%
  summarise(mean_cor = mean(partial_cor, na.rm = TRUE)) %>%
  ungroup()

# save all four correlations
cors = bind_rows(mutate(cors1, xval = 'mean vs. delta-variance'),
                 mutate(cors2, xval = 'variance vs. delta-variance'),
                 mutate(pcors1, xval = 'variance vs. delta-variance (partial)'),
                 mutate(pcors2, xval = 'mean vs. delta-variance (partial)')) %>%
  mutate(xval = fct_relevel(xval,
                            'variance vs. delta-variance (partial)',
                            'mean vs. delta-variance (partial)',
                            'variance vs. delta-variance',
                            'mean vs. delta-variance'))
saveRDS(cors, "data/analysis/mean_variance/correlations.rds")
