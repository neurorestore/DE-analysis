# Calculate the difference in variance between pseudobulks with biological and
# shuffled replicates across all 46 datasets.
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-write-delta-variance.R')
parser$add_argument('--input_file', type = 'character', required = T)
parser$add_argument('--output_dir', type = 'character', required = T)
args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(sparseMatrixStats)
source("R/functions/datasets.R")
source("R/functions/get_comparisons.R")

# set up output filepath
if (!dir.exists(args$output_dir))
  dir.create(args$output_dir, recursive = T)
dataset = args$input_file %>%
  basename() %>%
  gsub("\\.rds$", "", .)
output_filename = paste0(dataset, ".rds")
output_file = file.path(args$output_dir, output_filename)

# read input file and extract matrix/metadata
sc = readRDS(args$input_file)
expr = GetAssayData(sc, slot = 'counts')
meta = sc@meta.data
dataset = gsub("\\.rds$", "", basename(args$input_file))

# iterate through comparisons
vars = list()
comparisons = get_comparisons(dataset, expr, meta)
for (comparison_idx in seq_along(comparisons)) {
  comparison = comparisons[[comparison_idx]]
  comparison_name = names(comparisons)[comparison_idx]
  if (is.null(comparison_name))
    comparison_name = 1
  
  # get subset expression and metadata
  expr0 = comparison$expr
  meta0 = comparison$meta %>%
    set_rownames(colnames(expr0))
  
  # analyze each cell type in turn
  cell_types = unique(meta0$cell_type)
  for (cell_type_idx in seq_along(cell_types)) {
    cell_type = cell_types[cell_type_idx]
    message("  [", cell_type_idx, "/", length(cell_types),
            "] analyzing cell type: ", cell_type, " ...")
    
    # get cell-type-specific expression matrix
    keep = which(meta0$cell_type == cell_type)
    expr1 = expr0[, keep, drop = F]
    meta1 = meta0[keep, , drop = F]
    rownames(meta1) = colnames(expr1)
    
    # catch cell types without replicates or conditions
    if (n_distinct(meta1$label) < 2 | 
        n_distinct(meta1$replicate) < 3) {
      next
    }
    
    # shuffle replicates
    meta2 = meta1 %>%
      group_by(cell_type, label) %>%
      mutate(replicate = sample(replicate)) %>%
      ungroup() %>%
      set_rownames(colnames(expr1))
    
    # summarise to pseudobulk matrices
    metadatas = list('biological replicates' = meta1,
                     'shuffled replicates' = meta2)
    grid = tidyr::crossing(replicate_type = names(metadatas))
    
    tmp = data.frame()
    for (grid_idx in seq_len(nrow(grid))) {
      replicate_type = grid$replicate_type[grid_idx]
      model_matrix = metadatas[[replicate_type]] %>%
        ungroup() %>%
        mutate(label = as.character(label),
               replicate = as.character(replicate)) %>%
        model.matrix(~ 0 + replicate:label, data = .)
      mat_mm = expr1 %*% model_matrix
      
      # drop empty columns
      keep_samples = colSums(mat_mm) > 0
      mat_mm %<>% extract(, keep_samples) %>% as.matrix()
      
      # filter genes with 0 variance
      sds = rowSds(mat_mm)
      mat_mm %<>% extract(sds > 0, )
      
      # calculate CPM
      mat_mm %<>% edgeR::cpm()
      
      # calculate variances for each gene
      gene_vars = rowSds(mat_mm)
      
      # also calculate mean expression for each gene
      gene_means = rowMeans(mat_mm)
      
      # create output data frame
      df = data.frame(gene = rownames(mat_mm),
                      mean = gene_means,
                      variance = gene_vars,
                      replicate_type = replicate_type
      ) %>%
        # tag dataset, comparison, cell type
        mutate(dataset = dataset,
               comparison = comparison_name,
               cell_type = cell_type)
      tmp %<>% bind_rows(df)
    }
    
    # summarise within this cell type, at least
    summary = tmp %>%
      group_by(dataset, comparison, cell_type, gene) %>%
      filter(n() > 1) %>%
      mutate(cov = variance / mean,
             ratio = variance[replicate_type == 'shuffled replicates'] /
               variance[replicate_type == 'biological replicates'],
             delta = variance[replicate_type == 'shuffled replicates'] -
               variance[replicate_type == 'biological replicates'],
             delta_cov = cov[replicate_type == 'shuffled replicates'] -
               cov[replicate_type == 'biological replicates'],
             mean_var1 = mean(variance[replicate_type ==
                                         'shuffled replicates']),
             mean_var2 = mean(variance[replicate_type ==
                                         'biological replicates'])) %>%
      ungroup() %>%
      group_by(dataset, comparison, cell_type) %>%
      summarise(mean_ratio = mean(ratio, na.rm = TRUE),
                mean_delta = mean(delta, na.rm = TRUE),
                mean_delta_cov = mean(delta_cov, na.rm = TRUE),
                mean_var1 = mean(mean_var1, na.rm = TRUE),
                mean_var2 = mean(mean_var2, na.rm = TRUE)) %>%
      ungroup() %>%
      ## force everything to character
      map_dfc(as.character)
    
    # append
    vars %<>% bind_rows(summary)
  }
}

# save results
saveRDS(vars, output_file)
