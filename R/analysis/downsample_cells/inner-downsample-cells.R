# Run single-cell or pseudobulk DE analyses on downsampled datasets.
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-downsample-cells.R')
parser$add_argument('--input_file', type = 'character', required = T)
parser$add_argument('--output_dir', type = 'character', required = T)
parser$add_argument('--de_test', type = 'character', required = T)
parser$add_argument('--n_cells', type = 'double', required = T)
parser$add_argument('--sample_idx', type = 'integer', required = T)
args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
source("R/functions/run_DE.R")
source("R/functions/get_comparisons.R")

# set up output filepath
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = T)
dataset = args$input_file %>%
  basename() %>%
  gsub("\\.rds$", "", .)
output_filename = paste0(dataset,
                         "-de_test=", args$de_test,
                         "-n_cells=", args$n_cells,
                         "-sample_idx=", args$sample_idx,
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

  if (grepl("Hagai2018", dataset)) {
    # only do certain comparisons
    if (!comparison_name %in% c("lps4", "pic4"))
      next
  }

  message("[", comparison_idx, "/", length(comparisons), "] ",
          "analyzing comparison ", comparison_name, " ...")
  message("##############################")

  # get subset expression and metadata
  set.seed(args$sample_idx)
  meta0 = comparison$meta %>%
    # make sure rownames are correct
    set_rownames(colnames(comparison$expr)) %>%
    rownames_to_column(var = 'cell_barcode') %>%
    group_by(replicate) %>%
    mutate(cells = ceiling(args$n_cells * (n() / nrow(.)))) %>%
    sample_n(cells[1]) %>%
    ## maintaining the proportions, make sure n_cells is precise
    ungroup() %>%
    sample_n(args$n_cells) %>%
    set_rownames(.$cell_barcode)
  expr0 = comparison$expr %>% extract(, rownames(meta0))

  # make some checks
  if (grepl("Reyfman2020", dataset)) {
    cell_types = c("AT2", "Alveolar macrophages")
    meta0 %<>% filter(cell_type %in% cell_types) %>% set_rownames(.$cell_barcode)
    expr0 %<>% extract(, rownames(meta0))
  } else if (grepl("Angelidis2019", dataset)) {
        cell_types = c("Alveolar_macrophage", "Type_2_pneumocytes")
        meta0 %<>% filter(cell_type %in% cell_types) %>% set_rownames(.$cell_barcode)
        expr0 %<>% extract(, rownames(meta0))
  }

  # reconstruct the Seurat object
  sc_downsampled = CreateSeuratObject(expr0,
                                      min.cells = 1, min.features = 0,
                                      meta.data = meta0)

  # run DE analysis
  DE = run_DE(sc_downsampled, de_test = args$de_test)

  # append to list
  results[[comparison_name]] = DE
}

# stop if empty
if (length(results) == 0 | all(map_int(results, nrow) == 0))
  stop("couldn't get any results")

# save results
saveRDS(results, output_file)
