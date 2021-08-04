# Run gene set enrichment analysis (GSEA) on a set of DE results.
setwd("~/git/DE-analysis")
options(stringsAsFactors = FALSE)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-run-GSEA.R')
parser$add_argument('--input_file', type = 'character', required = TRUE)
parser$add_argument('--output_file', type = 'character', required = TRUE)
parser$add_argument('--n_permutations', type = 'integer', default = 1e6)
parser$add_argument('--min_size', type = 'integer', default = 10)
parser$add_argument('--max_size', type = 'integer', default = 1000)
args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(fgsea)
library(flavin)

# create output directory, if it does not exist
output_dir = dirname(args$output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = T)

# read input file
input = readRDS(args$input_file)

# read GO
species = ifelse(grepl("Angelidis|_mouse|_rat", args$input_file),
                 'mouse', 'human')
goa_file = paste0("data/GO/", 
                  fct_recode(species, 'mgi' = 'mouse', 'goa_human' = 'human'),
                  ".gaf.gz")
goa = read_gaf(goa_file)
ann = as_annotation_list(goa, 'DB_Object_Symbol', 'GO_ID')

# create results container
res = data.frame()

# iterate through comparisons
for (comparison_idx in seq_along(input)) {
  comparison = input[[comparison_idx]]
  comparison_name = names(input)[comparison_idx]
  if (is.null(comparison_name))
    comparison_name = 1
  
  if ("cell_type" %in% colnames(comparison)) {
    # iterate through cell types
    cell_types = unique(comparison$cell_type)
    ## keep only a subset of cell types to improve runtime
    keep = c("Naive",
             "Memory",
             "bone marrow derived mononuclear phagocytes",
             "Alveolar_macrophage",
             "Type_2_pneumocytes",
             "AT2",
             "Alveolar macrophages")
    cell_types %<>% intersect(keep)
    for (cell_type in cell_types) {
      message(".. analyzing cell type ", cell_type, " in comparison ",
              comparison_name, "...")
      DE = filter(comparison, cell_type == !!cell_type)
      
      # fix column names
      colnames(DE) %<>%
        fct_recode('p_val' = 'p.value',  ## DESeq2
                   'p_val' = 'pvalue',  ## DESeq2
                   'p_val' = 'p.value',  ## t/wilcox
                   'p_val' = 'P.Value',  ## limma
                   'p_val' = 'PValue'  , ## edgeR
                   'p_val_adj' = 'padj', ## DESeq2/t/wilcox
                   'p_val_adj' = 'adj.P.Val',      ## limma
                   'p_val_adj' = 'FDR',            ## edgeER
                   'avg_logFC' = 'log2FoldChange', ## DESEeq2
                   'avg_logFC' = 'logFC', ## limma/edgeR
                   'test_statistic' = 'stat', ## DESeq2
                   'test_statistic' = 'F', ## edgeR
                   'test_statistic' = 't', ## limma
                   'test_statistic' = 'LR', ## edgeR LRT
                   'test_statistic' = 'statistic' ## t
        ) %>%
        as.character()
      
      # extract ranks
      ranks = DE %>%
        drop_na(test_statistic) %$%
        setNames(abs(test_statistic), gene) %>% 
        sort(decreasing = TRUE)
      ## replace infinite values
      ranks[is.infinite(ranks)] = max(ranks[!is.infinite(ranks)])
      
      # run GSEA
      gsea = fgsea(pathways = ann,
                   stats = ranks,
                   nproc = 1,
                   nperm = args$n_permutations,
                   minSize = args$min_size,
                   maxSize = args$max_size) %>%
        dplyr::select(-leadingEdge) %>%
        # flag cell type and comparison
        mutate(cell_type = cell_type,
               comparison = comparison_name)
      
      # append to results
      res %<>% bind_rows(gsea)
    }
  } else {
    message(".. analyzing comparison ", comparison_name, "...")
    DE = comparison
    
    # fix column names
    colnames(DE) %<>%
      fct_recode('p_val' = 'p.value',  ## DESeq2
                 'p_val' = 'pvalue',  ## DESeq2
                 'p_val' = 'p.value',  ## t/wilcox
                 'p_val' = 'P.Value',  ## limma
                 'p_val' = 'PValue'  , ## edgeR
                 'p_val_adj' = 'padj', ## DESeq2/t/wilcox
                 'p_val_adj' = 'adj.P.Val',      ## limma
                 'p_val_adj' = 'FDR',            ## edgeER
                 'avg_logFC' = 'log2FoldChange', ## DESEeq2
                 'avg_logFC' = 'logFC', ## limma/edgeR
                 'test_statistic' = 'stat', ## DESeq2
                 'test_statistic' = 'F', ## edgeR
                 'test_statistic' = 't', ## limma
                 'test_statistic' = 'LR', ## edgeR LRT
                 'test_statistic' = 'statistic' ## t
      ) %>%
      as.character()
    
    # extract ranks
    ranks = DE %>%
      drop_na(test_statistic) %$%
      setNames(abs(test_statistic), gene) %>% 
      sort(decreasing = TRUE)
    
    # run GSEA
    gsea = fgsea(pathways = ann,
                 stats = ranks,
                 nproc = 1,
                 nperm = args$n_permutations,
                 minSize = args$min_size,
                 maxSize = args$max_size) %>%
      dplyr::select(-leadingEdge) %>%
      # flag cell type and comparison
      mutate(comparison = comparison_name)
    
    # append to results
    res %<>% bind_rows(gsea)
  }
}

# stop if empty
if (nrow(res) == 0)
  stop("couldn't get any results")

# save results
saveRDS(res, args$output_file)
