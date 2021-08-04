setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-run-concordance-terciles.R')
parser$add_argument('--input_sc', type = 'character', required = T)
parser$add_argument('--input_bulk', type = 'character', required = T)
parser$add_argument('--summary_file', type = 'character', required = T)
parser$add_argument('--output_dir', type = 'character', required = T)
parser$add_argument('--n_bins', type = 'integer', required = T)
args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(RRHO)
library(AUC)
source("R/functions/calculate_overlap.R")
source("R/analysis/bulk_concordance/write_grid.R")

# set up output filepath
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = T)

# load in files
sc = readRDS(args$input_sc)
bulk = readRDS(args$input_bulk)

# read expression summary
expr_summary = read.csv(args$summary_file)

# define output file
sc_name = gsub(".rds", "", basename(args$input_sc))
bulk_name = gsub(".rds", "", basename(args$input_bulk))
output_filename = paste0(sc_name, "|", bulk_name, '-n_bins=', args$n_bins,
                         ".rds")
output_file = file.path(args$output_dir, output_filename)

# get all combinations of single-cell/bulk
sc_idxs = names(sc)
bulk_idxs = names(bulk)
if (is.null(sc_idxs)) {
  sc_idxs = "1"
  names(sc) = '1'
}
if (is.null(bulk_idxs)) {
  bulk_idxs = "1"
  names(bulk) = '1'
}
comparisons = expand.grid(sc_comparison = sc_idxs, bulk_comparison = bulk_idxs,
                          stringsAsFactors = F)

# get rid of irrelevant comparisons from Cano-Gamez 2020
if (grepl("CanoGamez2020", sc_name)) {
  keep = map2_lgl(comparisons$sc_comparison,
                  comparisons$bulk_comparison,
                  ~ grepl(.x, .y))
  comparisons %<>% extract(keep, )
}

# analyze each comparison separately
results = data.frame()
for (comparison_idx in seq_len(nrow(comparisons))) {
  message("analyzing comparison ", comparison_idx, " of ", nrow(comparisons),
          " ...")
  
  # prepare data
  sc_comparison = comparisons$sc_comparison[comparison_idx]
  bulk_comparison = comparisons$bulk_comparison[comparison_idx]
  sc_sub = sc[[sc_comparison]]
  bulk_sub = bulk[[bulk_comparison]]
  comparison_label = paste0(sc_comparison, "|", bulk_comparison)
  
  # for Angelidis, filter to relevant cell types to prevent bugs
  if (grepl("Angelidis", sc_name)) {
    sc_sub %<>% filter(cell_type %in% c("Type_2_pneumocytes",
                                        "Alveolar_macrophage"))
  }
  # same for Reyfman
  if (grepl("Reyfman", sc_name)) {
    sc_sub %<>% filter(cell_type %in% c("AT2", "Alveolar macrophages"))
  }
  
  # run concordance over each quintile separately
  out = data.frame()
  cell_types = unique(sc_sub$cell_type)
  for (cell_type in cell_types) {
    # bin genes by expression level
    tested_genes = sc_sub %>%
      filter(cell_type == !!cell_type) %>%
      filter(gene %in% bulk_sub$gene) %>%
      pull(gene)
    bins = expr_summary %>%
      filter(gene %in% tested_genes) %>%
      filter(comparison == sc_comparison, cell_type == !!cell_type) %>%
      arrange(mean) %>%
      mutate(bin = cut(row_number() / n(),
                       breaks = seq(0, args$n_bins) / args$n_bins),
             bin = as.integer(bin)) %>%
      split(.$bin)
    
    # run over each bin
    for (bin in seq_len(args$n_bins)) {
      bin_genes = bins[[bin]]$gene
      sc_tmp = sc_sub %>%
        filter(cell_type == !!cell_type, gene %in% bin_genes)
      
      tmp = template %>%
        mutate(value = seq(nrow(template)) %>%
                 map( ~ {
                   print(template[., ])
                   method = template$method[.]
                   k = template$k[.]
                   cor_method = template$cor_method[.]
                   value = calculate_overlap(
                     bulk_de = bulk_sub,
                     sc_de = sc_tmp,
                     method = method,
                     k = k,
                     cor_method = cor_method
                   )
                 }) %>%
                 unlist()
        ) %>%
        mutate(cell_type = cell_type,
               bin = bin,
               sc_label = sc_comparison,
               bulk_label = bulk_comparison)
      out %<>% bind_rows(tmp)
    }
  }
  
  # append to the main results container
  results %<>% bind_rows(out)
}

# save results
saveRDS(results, output_file)
