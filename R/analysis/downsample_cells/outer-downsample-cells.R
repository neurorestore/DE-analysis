# Run single-cell or pseudobulk DE analyses on downsampled datasets.
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-downsample-cells.R')
parser$add_argument('--allocation', type = 'character')
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions/datasets.R")
source("R/functions/submit_job.R")
source("R/functions/detect_system.R")

# list input files
input_dir = file.path(base_dir, "rnaseq", "seurat")
input_files = file.path(input_dir, paste0(datasets, '.rds'))
inputs = data.frame(input_file = input_files)

# establish grid of analyses
opts = list(
  de_test = c(
    ## single-cell methods, implemented in Seurat
    "wilcox", "bimod", "t", "negbinom", "poisson", "LR", "MAST",
    ## pseudobulk methods
    "pseudobulk_DESeq2,test?LRT", "pseudobulk_DESeq2,test?Wald",
    "pseudobulk_limma,mode?voom", "pseudobulk_limma,mode?trend",
    "pseudobulk_edgeR,test?QLF", "pseudobulk_edgeR,test?LRT",
    ## mixed model methods, implemented in Seurat
    "mixed_lm", "mixed_nbinom", "mixed_poisson",
     ## slight adjustments to mixed model methods
    "mixed_lm,test?LRT", "mixed_nbinom,test?LRT", "mixed_poisson,test?LRT",
    "mixed_nbinom,offset?YES", "mixed_poisson,offset?YES",
    "mixed_nbinom,test?LRT,offset?YES", "mixed_poisson,test?LRT,offset?YES",
    ## pseudobulk methods with aggregation disabled
    "pseudobulk_DESeq2,test?LRT,replicate?cells",
    "pseudobulk_DESeq2,test?Wald,replicate?cells",
    "pseudobulk_limma,mode?voom,replicate?cells",
    "pseudobulk_limma,mode?trend,replicate?cells",
    "pseudobulk_edgeR,test?QLF,replicate?cells",
    "pseudobulk_edgeR,test?LRT,replicate?cells"
  ),
  n_cells = c(25, 50, 100, 200, 500, 1000),
  sample_idx = 0
)
grid = do.call(expand.grid, c(opts, stringsAsFactors = F))

# rep analysis grid over input files
grid %<>%
  dplyr::slice(rep(1:n(), each = nrow(inputs))) %>%
  mutate(input_file = rep(inputs$input_file, nrow(grid))) %>%
  left_join(inputs, by = 'input_file') %>%
  # reorder columns
  dplyr::select(input_file, de_test, n_cells, sample_idx) %>%
  # filter Angelidis when n_cells == 25 (won't run)
  filter(!grepl("Angelidis", input_file) | n_cells != 25)

# write the raw array
grid_file = "sh/analysis/downsample_cells/grids/downsample_cells.raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "analysis/downsample_cells/DE")

# check which parameters are already complete
overwrite = F
grid0 = grid
if (!overwrite) {
  grid0 = grid %>%
    mutate(idx = row_number(),
           output_file = file.path(output_dir,
                                   paste0(basename(input_file) %>%
                                            gsub("\\.rds$", "", .),
                                          '-de_test=', de_test,
                                          '-n_cells=', n_cells,
                                          '-sample_idx=', sample_idx,
                                          '.rds')),
           exists = file.exists(output_file)) %>%
    filter(!exists) %>%
    dplyr::select(-idx, -output_file, -exists)
}

# just do bulk datasets for now
grid0 %<>% filter(grepl("Hagai|Cano|Reyfman|Angelidis", input_file))

# write the grid that still needs to be run
write.table(grid0, "sh/analysis/downsample_cells/grids/downsample_cells.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
sh_dir = "~/git/DE-analysis/sh/analysis/downsample_cells"
script = file.path(sh_dir, "downsample_cells.sh")
submit_job(grid0, script, args$allocation, system)
