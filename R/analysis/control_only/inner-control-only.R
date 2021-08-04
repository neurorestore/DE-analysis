# Run single-cell or pseudobulk DE analyses on random splits of control
# samples only.
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-control-only.R')
parser$add_argument('--input_file', type = 'character', required = T)
parser$add_argument('--de_test', type = 'character', required = T)
parser$add_argument('--sample_idx', type = 'integer', required = T)
parser$add_argument('--shuffle_replicates', type = 'character', required = T)
parser$add_argument('--label', type = 'character', required = T)
parser$add_argument('--comparison', type = 'character', required = T)
parser$add_argument('--output_dir', type = 'character', required = T)
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
    '-de_test=', args$de_test,
    '-sample_idx=', args$sample_idx,
    '-shuffle_replicates=', args$shuffle_replicates,
    '-label=', args$label,
    '-comparison=', args$comparison,
  ".rds")
output_file = file.path(args$output_dir, output_filename)

# read input file and extract matrix/metadata
sc = readRDS(args$input_file)
expr = GetAssayData(sc, slot = 'counts')
meta = sc@meta.data

# get all combinations of conditions
comparisons = get_comparisons(dataset, expr, meta)
if (is.null(names(comparisons))) {
  names(comparisons) = '1'
}

# grab comparison of interest
comparison_name = args$comparison
comparison = comparisons[[comparison_name]]

# get subset expression and metadata
expr0 = comparison$expr
meta0 = comparison$meta

# grab the correct label
meta0 %<>%
  as.data.frame() %>%
  set_rownames(colnames(expr0)) %>%
  rownames_to_column(var = 'new_barcode') %>%
  filter(label == args$label) %>%
  set_rownames(.$new_barcode)
expr0 %<>% extract(, rownames(meta0))

# re-assign the groups
reps = unique(meta0$replicate)
n_reps = length(reps)
ctrl = sample(reps, n_reps/2)
meta0 %<>%
  mutate(label = ifelse(replicate %in% ctrl, 'ctrl', 'stim')) %>%
  set_rownames(.$new_barcode)

# check for replicate shuffling
if (args$shuffle_replicates == "YES") {
  meta0 %<>% group_by(cell_type, label) %>%
    mutate(replicate = sample(replicate)) %>%
    set_rownames(.$new_barcode)
}

# reconstruct the Seurat object
sc0 = CreateSeuratObject(expr0, min.cells = 1, min.features = 0,
                         meta.data = meta0)

# run DE analysis
DE = run_DE(sc0, de_test = args$de_test)

# stop if empty
if (nrow(DE) == 0)
  stop("couldn't get any results")

# save results
saveRDS(DE, output_file)
