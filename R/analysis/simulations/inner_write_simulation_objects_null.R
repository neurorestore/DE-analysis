# Generate the complete set of simulated scRNA-seq datasets for the experiment
# of null differential expression
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner_write_simulation_objects_null.R')
parser$add_argument('--n_cells', type = 'integer', required = T)
parser$add_argument('--de_prob', type = 'double', required = T)
parser$add_argument('--de_loc', type = 'double', required = T)
parser$add_argument('--n_reps', type = 'integer', required = T)
parser$add_argument('--sample_idx', type = 'integer', required = T)
parser$add_argument('--output_dir', type = 'character', required = T)
args = parser$parse_args()
print(args)

library(Seurat)
library(splatterBatch)
library(scater)
library(tidyverse)
library(magrittr)
library(Matrix)

source("R/functions/detect_system.R")

# check the output directory
if (!dir.exists(args$output_dir)) {
  dir.create(args$output_dir, recursive = T)
}

# define output file
output_filename = paste0("GSE96583",
                         "-n_cells=", args$n_cells,
                         "-de_prob=", args$de_prob,
                         "-de_loc=", args$de_loc,
                         "-n_reps=", args$n_reps,
                         "-sample_idx=", args$sample_idx,
                         ".rds")
output_file = file.path(args$output_dir, output_filename)

# get parameters defined by Kang et al. IFN dataset
params = readRDS(file.path(base_dir, "analysis/simulations/parameters",
                           "parameters_GSE96583.rds"))

# calculate group probabilities
group_probs = 1 / args$n_reps
# assign groups
unst = sample(seq(args$n_reps), args$n_reps/2)

# generate simulated cells
sim = splatterBatch::splatSimulateGroups(
  params = params,
  seed = args$sample_idx,
  batchCells = args$n_cells,
  de.prob = args$de_prob,
  de.facLoc = args$de_loc,
  group.prob = rep(group_probs, args$n_reps), verbose = F
) %>% logNormCounts() %>% as.Seurat()

# adjust metadata for default input to Seurat
sim@meta.data %<>%
  dplyr::mutate(cell_type = paste0('cell_', 1)) %>%
  dplyr::rename(label = Group) %>%
  mutate(replicate = gsub("Group", "Replicate ", label)) %>%
  mutate(label = as.numeric(gsub("Group", "", label))) %>%
  mutate(label = ifelse(label %in% unst, 'unst', 'stim')) %>%
  set_rownames(colnames(GetAssayData(sim)))

# save
saveRDS(sim, output_file)
