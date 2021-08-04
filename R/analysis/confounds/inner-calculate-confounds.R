# Calculate possible confounding factors to DE identified in the Augur paper:
# - read depth (mean and total counts per cell type),
# - 'gene depth' (number of genes detected per cell type), and
# - 'cell depth' (number of cells sequenced per type)
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-calculate-confounds.R')
parser$add_argument('--input_file', type = 'character', required = T)
parser$add_argument('--output_dir', type = 'character', required = T)
args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
source("R/functions/get_comparisons.R")

# set up output filepath
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = T)
dataset = args$input_file %>%
  basename() %>% 
  gsub("\\.rds$", "", .)
output_filename = paste0(dataset, ".txt")
output_file = file.path(args$output_dir, output_filename)

# read input file and extract matrix/metadata
sc = readRDS(args$input_file)
expr = GetAssayData(sc, slot = 'counts')
meta = sc@meta.data
dataset = gsub("\\.rds$", "", basename(args$input_file))

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
  meta0 = comparison$meta %>%
    set_rownames(colnames(expr0))
  
  # reconstruct the Seurat object
  sc0 = CreateSeuratObject(expr0, min.cells = 1, min.features = 0,
                           meta.data = meta0)
  
  # analyze each cell type in turn
  cell_types = unique(meta0$cell_type) 
  for (cell_type_idx in seq_along(cell_types)) {
    cell_type = cell_types[cell_type_idx]
    message("  [", cell_type_idx, "/", length(cell_types), 
            "] analyzing cell type: ", cell_type, " ...")
    
    # number of cells
    keep = which(meta0$cell_type == cell_type)
    n_cells = length(keep)
    
    # read depth per cell
    expr1 = expr0[, keep, drop = F]
    reads = colSums(expr1)
    read_depth_mean = mean(reads)
    read_depth_sum = sum(reads)
    
    # genes detected per cell
    n_genes = colSums(expr1 > 0)
    n_genes_mean = mean(n_genes)

    # append to results
    results %<>% rbind(
      data.frame(dataset = dataset,
                 comparison = comparison_name,
                 cell_type = cell_type,
                 outcome = c("# of cells",
                             "read depth (mean)",
                             "read depth (sum)",
                             "# of genes (mean)"),
                 value = c(n_cells,
                           read_depth_mean, 
                           read_depth_sum,
                           n_genes_mean)))
  }
}

# write
write.csv(results, output_file, row.names = F)
system(paste("gzip --force", output_file))
