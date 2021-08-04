# Run single-cell or pseudobulk DE analyses on the Hagai et al. dataset with 
# ERCC spike-ins.
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-run-DE.R')
parser$add_argument('--allocation', type = 'character')
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions/submit_job.R")
source("R/functions/detect_system.R")

# list input files
input_dir = file.path(base_dir, "rnaseq", "seurat")
input_files = file.path(input_dir, "Hagai2018_plate.rds")
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
    ## mixed model, implemented in Seurat
    "mixed_lm",
    ## pseudobulk methods run without aggregation
    "pseudobulk_DESeq2,test?LRT,replicate?cells",
    "pseudobulk_DESeq2,test?Wald,replicate?cells",
    "pseudobulk_limma,mode?voom,replicate?cells",
    "pseudobulk_limma,mode?trend,replicate?cells",
    "pseudobulk_edgeR,test?QLF,replicate?cells",
    "pseudobulk_edgeR,test?LRT,replicate?cells"
  ),
  shuffle_replicates = c("NO", "YES")
)
grid = do.call(expand.grid, c(opts, stringsAsFactors = F))

# rep analysis grid over input files
grid %<>%
  dplyr::slice(rep(1:n(), each = nrow(inputs))) %>%
  mutate(input_file = rep(inputs$input_file, nrow(grid))) %>%
  left_join(inputs, by = 'input_file') %>%
  # reorder columns
  dplyr::select(input_file, everything())

# write the raw array
grid_file = "sh/analysis/run_spike_in_DE/grids/run_DE.raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "analysis/run_spike_in_DE")

# check which parameters are already complete
overwrite = F
grid0 = grid
if (!overwrite) {
  grid0 = grid %>%
    mutate(output_filename = paste0(basename(input_file) %>%
                                      gsub("\\.rds$", "", .),
                                    '-de_test=', de_test,
                                    '-shuffle_replicates=', shuffle_replicates,
                                    '.rds'),
           output_file = file.path(output_dir, output_filename),
           exists = file.exists(output_file)) %>%
    filter(!exists) %>%
    dplyr::select(-output_file, -output_filename, -exists)
}

# subset grid, if needed
if (nrow(grid0) >= 10000) {
  grid0 %<>% dplyr::slice(1:9900) ## allow for some other running jobs or sh
}

# write the grid that still needs to be run
write.table(grid0, "sh/analysis/run_spike_in_DE/grids/run_DE.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
sh_dir = "~/git/DE-analysis/sh/analysis/run_spike_in_DE"
script = file.path(sh_dir, "run_DE.sh")
submit_job(grid0, script, args$allocation, system)
