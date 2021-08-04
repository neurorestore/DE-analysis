# Run bulk DE analyses on all cell types in a dataset.
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-run-DE.R')
parser$add_argument('--input_file', type = 'character', required = T)
parser$add_argument('--output_dir', type = 'character', required = T)
parser$add_argument('--de_test', type = 'character', required = T)
args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
source("R/functions/get_bulk_comparisons.R")
source("R/functions/run_DE.R")

# set up output filepath
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = T)
dataset = args$input_file %>%
  basename() %>%
  gsub("\\.rds$", "", .)
dataset_label = args$input_file %>%
  basename() %>%
  gsub("_.*|.rds", "", .)
output_filename = paste0(dataset, "-de_test=", args$de_test, ".rds")
output_file = file.path(args$output_dir, output_filename)

# read input file and extract matrix/metadata
sc = readRDS(args$input_file)
expr = sc$assay
meta = sc$meta

# get all combinations of conditions
results = list()
comparisons = get_bulk_comparisons(dataset_label, expr, meta)
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
  meta0 = comparison$meta %>%
    set_rownames(colnames(expr0))

  # run DE analysis
  if (grepl("proteomics|microarray", args$input_file)) {
    DE = bulk_DE(expr0, targets = meta0, de_test = args$de_test, used_voom = F)
  } else {
    DE = bulk_DE(expr0, targets = meta0, de_test = args$de_test)
  }
  # append to list
  results[[comparison_name]] = DE
}

# stop if empty
if (length(results) == 0 | all(map_int(results, nrow) == 0))
  stop("couldn't get any results")

# save results
saveRDS(results, output_file)
