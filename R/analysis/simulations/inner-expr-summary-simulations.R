setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-expr-summary.R')
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
output_filename = basename(args$input_file)
output_file = file.path(args$output_dir, output_filename)

# read input file and extract matrix/metadata
sc = readRDS(args$input_file)
expr = GetAssayData(sc, slot = 'counts')
meta = sc@meta.data

# calculate statistics
genes = rownames(expr)
means = Matrix::rowMeans(expr)
sds = sparseMatrixStats::rowSds(expr)
covs = sds / means
pct_zeros = Matrix::rowSums(expr == 0) / ncol(expr)

# calculate logFC as defined in Seurat
logFC = tryCatch({
  sc0 = CreateSeuratObject(expr, meta = meta) %>%
    NormalizeData()
  Idents(sc0) = sc0$label
  mat = GetAssayData(sc0, slot = 'data')
  levels = levels(meta$label)
  if (is.null(levels)) {
    levels = unique(meta$label)
  }
  cells1 = WhichCells(sc0, idents = levels[1])
  cells2 = WhichCells(sc0, idents = levels[2])
  data1 = log(rowMeans(mat[, cells1, drop = F] + 1))
  data2 = log(rowMeans(mat[, cells2, drop = F] + 1))
  out = data2 - data1 # backwards from Seurat (i.e., the proper way)
}, error = function(e) { return(NA_real_) })

# calculate pseudobulk variance
pseudobulk_variance = tryCatch({
  meta2 = meta %>%
    mutate(label = as.character(label),
           replicate = as.character(replicate))
  mm = model.matrix(~ 0 + replicate, data = meta2)
  mat_mm = expr %*% mm
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
  meta2 = meta %>%
    mutate(label = as.character(label),
           replicate = as.character(replicate)) %>%
    group_by(cell_type, label) %>%
    mutate(replicate = sample(replicate))
  mm = model.matrix(~ 0 + replicate, data = meta2)
  mat_mm = expr %*% mm
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
  # drop genes with zero expression
  filter(mean > 0)

# write
saveRDS(df, output_file)
