# Calculate outcomes for single-cell or pseudobulk DE analyses on
# downsampled datasets.
setwd("~/git/DE-analysis")
options(stringsAsFactors = F)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'outer-downsample-cells-outcomes.R')
parser$add_argument('--allocation', type = 'character')
args = parser$parse_args()

library(tidyverse)
library(magrittr)
source("R/functions/datasets.R")
source("R/functions/submit_job.R")
source("R/functions/detect_system.R")

# list bulk input files
bulk_files = list.files(file.path(base_dir, "analysis/run_bulk_DE"))
bulk_inputs = data.frame(bulk_file = bulk_files) %>%
  mutate(label = gsub("-.*", "", bulk_file)) %>%
  # manual fix for the Hagai datasets
  mutate(label = ifelse(label == 'Hagai2018', gsub("-.*", "", bulk_file),
                        label)) %>%
  # manually match a few sc datasets to their bulk data
  mutate(label = fct_recode(label,
                            "Reyfman2020" = "Reyfman2020_alvmac",
                            "Reyfman2020" = "Reyfman2020_AT2",
                            "Angelidis2019" = "Angelidis2019_facsepi",
                            "Angelidis2019" = "Angelidis2019_facsmac",
                            "CanoGamez2020" = "CanoGamez2020:proteomics")) %>%
  # restore the entire filepath
  mutate(bulk_file = file.path(base_dir, 'analysis/run_bulk_DE', bulk_file))

# list datasets
inputs = data.frame(dataset = datasets)

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
    "pseudobulk_edgeR,test?LRT,replicate?cells"),
  n_cells = c(25, 50, 100, 200, 500, 1000),
  sample_idx = 0
)
grid = do.call(expand.grid, c(opts, stringsAsFactors = F))

# rep analysis grid over input files
grid %<>%
  dplyr::slice(rep(1:n(), each = nrow(inputs))) %>%
  mutate(dataset = rep(inputs$dataset, nrow(grid))) %>%
  left_join(inputs, by = 'dataset') %>%
  # reorder columns
  dplyr::select(dataset, de_test, n_cells, sample_idx) %>%
  mutate(label = dataset) %>%
  # add in bulk file
  left_join(bulk_inputs) %>%
  # filter Angelidis when n_cells == 25 (won't run)
  filter(!grepl("Angelidis", dataset) | n_cells != 25)

# now, reorganize the grid to map query/target file pairs
query_dir = file.path(base_dir, "analysis/downsample_cells/DE")
grid %<>%
  mutate(input_sc = file.path(query_dir,
                                paste0(dataset,
                                       '-de_test=', de_test,
                                       '-n_cells=', n_cells,
                                       '-sample_idx=', sample_idx,
                                       '.rds'))) %>%
  dplyr::rename(input_bulk = bulk_file) %>%
  dplyr::select(input_sc, input_bulk)

# write the raw array
grid_file = "sh/analysis/downsample_cells/grids/downsample_cells.raw.txt"
grid_dir = dirname(grid_file)
if (!dir.exists(grid_dir))
  dir.create(grid_dir, recursive = T)
write.table(grid, grid_file, quote = F, row.names = F, sep = "\t")

# define output directory where results are stored
output_dir = file.path(base_dir, "analysis/downsample_cells/concordance")

# check which parameters are already complete
overwrite = F
grid0 = grid
if (!overwrite) {
  grid0 = grid %>%
    mutate(idx = row_number(),
           output_filename = paste0(basename(input_sc) %>%
                                             gsub("\\.rds$", "", .),
                                           "|",
                                           basename(input_bulk) %>%
                                             gsub("\\.rds$", "", .),
                                           '.rds'),
                  output_file = file.path(output_dir, output_filename),
                  exists = file.exists(output_file)) %>%
    filter(!exists) %>%
    dplyr::select(-idx, -output_file, -exists)
}

# just do hagai and Canogamez for now
grid0 %<>% filter(grepl("Hagai|Cano|Reyfman|Angelidis", input_sc))

# write the grid that still needs to be run
write.table(grid0,
            "sh/analysis/downsample_cells/grids/downsample_cells_outcomes.txt",
            quote = F, row.names = F, sep = "\t")

# finally, run the job on whatever system we're on
sh_dir = "~/git/DE-analysis/sh/analysis/downsample_cells"
script = file.path(sh_dir, "downsample_cells_outcomes.sh")
submit_job(grid0, script, args$allocation, system)
