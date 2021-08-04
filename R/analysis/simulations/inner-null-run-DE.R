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

# check for replicate shuffling
if (args$shuffle_replicates == "YES") {
  sc@meta.data %<>%
    group_by(label) %>%
    mutate(replicate = sample(replicate)) %>%
    ungroup() %>%
    set_rownames(colnames(sc))
}

# run DE analysis
DE = run_DE(sc, de_test = args$de_test)

# stop if empty
if (nrow(DE) == 0)
  stop("couldn't get any results")

# save results
saveRDS(DE, output_file)
