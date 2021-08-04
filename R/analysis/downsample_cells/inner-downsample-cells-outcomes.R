# Get concordance between scRNA-seq/pseudobulk DE and bulk DE
# in downsampled data
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-downsample-cells-outcomes.R')
parser$add_argument('--input_sc', type = 'character', required = T)
parser$add_argument('--input_bulk', type = 'character', required = T)
parser$add_argument('--output_dir', type = 'character', required = T)
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

# define output file
sc_name = gsub(".rds", "", basename(args$input_sc))
bulk_name = gsub(".rds", "", basename(args$input_bulk))
output_filename = paste0(sc_name, "|", bulk_name, ".rds")
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
comparisons = expand.grid(sc_idxs, bulk_idxs, stringsAsFactors = F)

# get rid of irrelevant comparisons from Cano-Gamez 2020
if (grepl("CanoGamez2020", sc_name)) {
  keep = map2_lgl(comparisons$Var1, comparisons$Var2, ~ grepl(.x, .y))
  comparisons %<>% extract(keep, )
}

results = c()
for (comparison_idx in 1:nrow(comparisons)) {
  message("analyzing comparison ", comparison_idx, " of ", nrow(comparisons),
          " ...")

  # prepare data
  sc_sub = sc[[comparisons[comparison_idx,]$Var1]]
  bulk_sub = bulk[[comparisons[comparison_idx,]$Var2]] %>% 
    ## fix for Reyfman
    ungroup()
  sc_label = comparisons[comparison_idx,]$Var1
  bulk_label = comparisons[comparison_idx,]$Var2
  comparison_label = paste0(sc_label, "|", bulk_label)
  
  # for Angelidis, filter to relevant cell types to prevent bugs
  if (grepl("Angelidis", sc_name)) {
    sc_sub %<>% filter(cell_type %in% c("Type_2_pneumocytes",
                                        "Alveolar_macrophage"))
  }
  # same for Reyfman
  if (grepl("Reyfman", sc_name)) {
    sc_sub %<>% filter(cell_type %in% c("AT2", "Alveolar macrophages"))
  }
  
  # calculate concordance metrics for this comparison
  out = sc_sub %>%
    split(.$cell_type) %>%
    map( ~ {
      print(.$cell_type[1])
      sc_tmp = .
      tmp = template %>%
        mutate(value = seq(nrow(template)) %>%
                 map( ~ {
                   print(template[., ])
                   method = template[.,]$method
                   k = template[.,]$k
                   cor_method = template[.,]$cor_method
                   value = calculate_overlap(
                     bulk_de = bulk_sub,
                     sc_de = sc_tmp,
                     method = method,
                     k = k,
                     cor_method = cor_method)
                 }) %>%
                 unlist()
        )
    }) %>%
    bind_rows(.id = 'cell_type') %>%
    mutate(
      sc_label = sc_label,
      bulk_label = bulk_label
    )
  # bind to main results container
  results %<>% bind_rows(out)
}

# save results
saveRDS(results, output_file)
