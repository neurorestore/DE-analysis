# Run single-cell or pseudobulk DE analyses on all cell types in a dataset.
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-run-DE.R')
parser$add_argument('--input_file', type = 'character', required = T)
parser$add_argument('--shuffle_replicates', type = 'character', required = T)
parser$add_argument('--output_dir', type = 'character', required = T)
parser$add_argument('--de_test', type = 'character', required = T)
args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(peakRAM)
library(future)
source("R/functions/get_comparisons.R")
source("R/functions/run_DE.R")

# set up output filepath
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = T)
dataset = args$input_file %>%
  basename() %>%
  gsub("\\.rds$", "", .)
output_filename = paste0(dataset,
                         "-de_test=", args$de_test,
                         "-shuffle_replicates=", args$shuffle_replicates,
                         ".rds")
output_file = file.path(args$output_dir, output_filename)

# read input file and extract matrix/metadata
sc = readRDS(args$input_file)
expr = GetAssayData(sc, slot = 'counts')
meta = sc@meta.data

# get all combinations of conditions
results = list()
comparisons = get_comparisons(dataset, expr, meta)
for (comparison_idx in seq_along(comparisons)) {
  comparison = comparisons[[comparison_idx]]
  comparison_name = names(comparisons)[comparison_idx]
  if (is.null(comparison_name))
    comparison_name = 1
  
  message("[", comparison_idx, "/", length(comparisons), "] ",
          "analyzing comparison ", comparison_name, " ...")
  message("##############################")
  
  # get subset expression and metadata
  expr0 = comparison$expr
  meta0 = comparison$meta
  
  # check for replicate shuffling
  if (args$shuffle_replicates == "YES") {
    meta0 %<>% 
      group_by(cell_type, label) %>%
      mutate(replicate = sample(replicate))
  }
  
  # fix rownames
  meta0 %<>% set_rownames(colnames(expr0))
  
  # reconstruct the Seurat object
  sc0 = CreateSeuratObject(expr0, min.cells = 1, min.features = 0,
                           meta.data = meta0)
  
  # run DE analysis
  DE = run_DE(sc0, de_test = args$de_test)
  
  # append to list
  results[[comparison_name]] = DE
}

# stop if empty
if (length(results) == 0 | all(map_int(results, nrow) == 0))
  stop("couldn't get any results")

# save results
saveRDS(results, output_file)
