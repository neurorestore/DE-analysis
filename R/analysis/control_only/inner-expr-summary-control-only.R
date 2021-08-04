setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-expr-summary.R')
parser$add_argument('--input_file', type = 'character', required = T)
parser$add_argument('--label', type = 'character', required = T)
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
output_filename = paste0(dataset, "-label=", args$label, ".csv")
output_file = file.path(args$output_dir, output_filename)

# read input file and extract matrix/metadata
sc = readRDS(args$input_file)
expr = GetAssayData(sc, slot = 'counts')
meta = sc@meta.data
dataset = gsub("\\.rds$", "", basename(args$input_file))

# filter to the condition of interest
keep = which(meta$label == args$label) 
expr0 = expr[, keep, drop = F]
meta0 = meta[keep, , drop = F]

# analyze each cell type in turn
results = data.frame()
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
  
  # calculate statistics
  genes = rownames(expr1)
  means = Matrix::rowMeans(expr1)
  sds = sparseMatrixStats::rowSds(expr1)
  covs = sds / means
  pct_zeros = Matrix::rowSums(expr1 == 0) / ncol(expr1)
  
  # calculate logFC as defined in Seurat
  logFC = tryCatch({
    sc0 = CreateSeuratObject(expr1, meta = meta1) %>%
      NormalizeData()
    Idents(sc0) = sc0$label
    mat = GetAssayData(sc0, slot = 'data')
    levels = levels(meta1$label)
    if (is.null(levels)) {
      levels = unique(meta1$label)
    }
    cells1 = WhichCells(sc0, idents = levels[1])
    cells2 = WhichCells(sc0, idents = levels[2])
    data1 = log(rowMeans(mat[, cells1, drop = F] + 1))
    data2 = log(rowMeans(mat[, cells2, drop = F] + 1))
    out = data2 - data1 # backwards from Seurat (i.e., the proper way)
  }, error = function(e) { return(NA_real_) })
  
  # calculate pseudobulk variance
  pseudobulk_variance = tryCatch({
    meta2 = meta1 %>%
      mutate(label = as.character(label),
             replicate = as.character(replicate))
    mm = model.matrix(~ 0 + replicate, data = meta2)
    mat_mm = expr1 %*% mm
    # drop empty columns
    keep_samples = colSums(mat_mm) > 0
    mat_mm %<>% extract(, keep_samples) %>% as.matrix()
    # normalize
    mat_mm %<>% edgeR::cpm()
    # grab the variance for each gene
    vars = sparseMatrixStats::rowSds(mat_mm)
    vars %<>% setNames(rownames(mat_mm))
    vars
  }, error = function(e) { return(NA_real_) })
  
  # calculate shuffled pseudobulk variance
  shuffled_variance = tryCatch({
    meta2 = meta1 %>%
      mutate(label = as.character(label),
             replicate = as.character(replicate)) %>%
      group_by(cell_type, label) %>%
      mutate(replicate = sample(replicate))
    mm = model.matrix(~ 0 + replicate, data = meta2)
    mat_mm = expr1 %*% mm
    # drop empty columns
    keep_samples = colSums(mat_mm) > 0
    mat_mm %<>% extract(, keep_samples) %>% as.matrix()
    # normalize
    mat_mm %<>% edgeR::cpm()
    # grab the variance for each gene
    vars = sparseMatrixStats::rowSds(mat_mm)
    vars %<>% setNames(rownames(mat_mm))
    vars
  }, error = function(e) { return(NA_real_) })
  
  # calculate the ratio of real to shuffled variance
  ratio = pseudobulk_variance / shuffled_variance
  
  # convert to data frame
  df = data.frame(gene = genes, mean = means, sd = sds, cov = covs,
                  pct_zero = pct_zeros, logFC = logFC,
                  pseudobulk_variance = pseudobulk_variance,
                  shuffled_variance = shuffled_variance,
                  pseudobulk_ratio = ratio) %>%
    mutate(dataset = dataset,
           label = args$label,
           cell_type = cell_type)
  
  # append to results
  results %<>% bind_rows(df)
}

# rearrange columns
results %<>% dplyr::select(dataset, label, cell_type, everything())

# write
write.csv(results, output_file, row.names = F)
system(paste("gzip --force", output_file))
