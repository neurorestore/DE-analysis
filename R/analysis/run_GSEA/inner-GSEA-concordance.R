# Calculate the concordance GSEA results from matching single-cell and bulk DE.
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-GSEA-concordance.R')
parser$add_argument('--label', type = 'character', required = TRUE)
parser$add_argument('--input_sc', type = 'character', required = TRUE)
parser$add_argument('--input_bulk', type = 'character', required = TRUE)
parser$add_argument('--output_file', type = 'character', required = TRUE)
args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(AUC)
source("R/functions/calculate_overlap.R")
source("R/analysis/bulk_concordance/write_grid.R")

# set up output filepath
output_dir = dirname(args$output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)

# load in files
sc = readRDS(args$input_sc)
bulk = readRDS(args$input_bulk)

# iterate through single-cell comparisons
res = data.frame()
label = args$label
for (sc_comparison in unique(sc$comparison)) {
  sc_sub = filter(sc, comparison == sc_comparison)
  
  # iterate through cell types in the single-cell data
  cell_types = unique(sc_sub$cell_type)
  for (cell_type in cell_types) {
    message(".. analyzing cell type ", cell_type, " in comparison ",
            sc_comparison, "...")
    input1 = filter(sc_sub, cell_type == !!cell_type)
    
    # now, get the matching bulk data
    if (grepl("Hagai2018", label)) {
      bulk_comparison = toupper(sc_comparison)
      input2 = filter(bulk, comparison == bulk_comparison)
    } else if (label == "CanoGamez2020") {
      bulk_comparison = paste0('Resting|', sc_comparison, '|', cell_type, '|5d')
      input2 = filter(bulk, comparison == bulk_comparison)
    } else if (grepl("Reyfman2020|Angelidis2019", label)) {
      bulk_comparison = 1
      input2 = bulk
    } else {
      stop("not sure what to do with label: ", label)
    }
    
    # fix columns
    input1 %<>% dplyr::rename(p_val = pval, p_val_adj = padj, gene = pathway,
                              test_statistic = nMoreExtreme) %>%
      mutate(avg_logFC = 1) ## need to set the sign
    input2 %<>% dplyr::rename(p_val = pval, p_val_adj = padj, gene = pathway,
                              test_statistic = nMoreExtreme) %>%
      mutate(avg_logFC = 1)
    
    # run the GSEA results through our generic concordance function
    concordance = template %>%
      mutate(value = pmap_dbl(., function(...) {
        template = tibble(...)
        print(template)
        value = calculate_overlap(
          bulk_de = input2,
          sc_de = input1,
          method = template$method,
          k = template$k,
          cor_method = template$cor_method
        )
        return(value)
      })) %>%
      # flag comparisons and cell type
      mutate(sc_comparison = sc_comparison,
             bulk_comparison = bulk_comparison,
             cell_type = cell_type)
    
    # append to results
    res %<>% rbind(concordance)
  }
}

# save results
saveRDS(res, args$output_file)
