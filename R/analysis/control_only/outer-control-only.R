# Run single-cell or pseudobulk DE analyses on random splits of control
# samples only.
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-control-only.R')
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
    # pseudobulk methods
    "pseudobulk_DESeq2,test?LRT", "pseudobulk_DESeq2,test?Wald",
    "pseudobulk_limma,mode?voom", "pseudobulk_limma,mode?trend",
    "pseudobulk_edgeR,test?QLF",
    "pseudobulk_edgeR,test?LRT",
    ## mixed model methods, implemented in Seurat
    "mixed_lm"
  ),
  sample_idx = 1,
  shuffle_replicates = c("NO", "YES")
)
grid = do.call(expand.grid, c(opts, stringsAsFactors = F)) %>%
  mutate(type = ifelse(grepl("pseudo|mixed", de_test), 'rep', 'single')) %>%
  # only do shuffle replicates in pseudobulk or mixed model
  filter(type != 'single' | shuffle_replicates != 'YES') %>%
  dplyr::select(-type)

# rep analysis grid over input files
grid %<>%
  dplyr::slice(rep(1:n(), each = nrow(inputs))) %>%
  mutate(input_file = rep(inputs$input_file, nrow(grid))) %>%
  left_join(inputs, by = 'input_file') %>%
  # reorder columns
  dplyr::select(input_file, de_test, sample_idx, shuffle_replicates)

# load in number of replicates for each dataset
reps = readRDS("data/analysis/confounds/replicates.rds")

# grab the conditions where we have enough replicates
keep = reps %>%
  group_by(dataset, label, comparison) %>%
  summarise(n_reps = n()) %>%
  ungroup %>%
  filter(n_reps >= 6)

# add this into the grid
grid %<>%
    mutate(dataset = gsub(".rds$", "", basename(input_file))) %>%
    left_join(keep) %>%
    # drop the ones we aren't keeping
    drop_na() %>%
    dplyr::select(-dataset, -n_reps)

# write the raw array
grid_file = "sh/analysis/control_only/grids/control_only.raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "analysis/control_only")

# check which parameters are already complete
overwrite = F
grid0 = grid
if (!overwrite) {
  grid0 = grid %>%
    mutate(output_file = file.path(output_dir, paste0(basename(input_file) %>%
                                          gsub("\\.rds$", "", .),
                                        '-de_test=', de_test,
                                        '-sample_idx=', sample_idx,
                                        '-shuffle_replicates=', shuffle_replicates,
                                        '-label=', label,
                                        '-comparison=', comparison,
                                        '.rds')),
           exists = file.exists(output_file)) %>%
    filter(!exists) %>%
    dplyr::select(-output_file, -exists)
}

# write the grid that still needs to be run
write.table(grid0, "sh/analysis/control_only/grids/control_only.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
sh_dir = "~/git/DE-analysis/sh/analysis/control_only"
script = file.path(sh_dir, "control_only.sh")
submit_job(grid0, script, args$allocation, system)
