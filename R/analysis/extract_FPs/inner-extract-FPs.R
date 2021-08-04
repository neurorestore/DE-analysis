# Extract summary statistics for the top-ranking false-positives and
# false-negatives from each DE method.
setwd("~/git/DE-analysis")
options(stringsAsFactors = FALSE)
library(argparse)

# parse arguments
parser = ArgumentParser(prog = 'inner-extract-FPs.R')
parser$add_argument('--label', type = 'character', required = TRUE)
parser$add_argument('--sc_file', type = 'character', required = TRUE)
parser$add_argument('--bulk_file', type = 'character', required = TRUE)
parser$add_argument('--summary_file', type = 'character', required = TRUE)
parser$add_argument('--output_file', type = 'character', required = TRUE)
args = parser$parse_args()

library(tidyverse)
library(magrittr)

# read single-cell and bulk DE
sc_de = readRDS(args$sc_file)
bulk_de = readRDS(args$bulk_file)

# read expr summary
expr_summary = read.csv(args$summary_file)

# set up output containers
FP = FN = data.frame()

# iterate through single-cell comparisons
label = args$label
for (comparison_idx in seq_along(sc_de)) {
  sc_sub = sc_de[[comparison_idx]]
  sc_comparison = names(sc_de)[comparison_idx]
  if (is.null(sc_comparison))
    sc_comparison = 1
  
  # filter comparisons
  if (grepl("Hagai2018", label) & !sc_comparison %in% c("lps4", "pic4")) {
    message(".. skipping comparison ", sc_comparison, "...")
    next
  }

  # iterate through cell types in the single-cell data
  cell_types = unique(sc_sub$cell_type)
  ## for Reyfman/Angelidis, only do select cell types
  if (grepl("Reyfman2020", label)) {
    cell_types = ifelse(grepl("AT2", label), "AT2", "Alveolar macrophages")
  } else if (grepl("Angelidis2019", label)) {
    cell_types = ifelse(grepl("alvmac", label), "Alveolar_macrophage",
                        "Type_2_pneumocytes")
  }
  for (cell_type in cell_types) {
    message(".. analyzing cell type ", cell_type, " in comparison ",
            sc_comparison, "...")
    sc = filter(sc_sub, cell_type == !!cell_type)
    
    # get the relevant bulk comparisons
    if (grepl("Hagai2018", label)) {
      bulk_comparison = toupper(sc_comparison)
      bulk = bulk_de[[bulk_comparison]]
    } else if (label == "CanoGamez2020") {
      bulk_comparison = paste0('Resting|', sc_comparison, '|', cell_type, '|5d')
      bulk = bulk_de[[bulk_comparison]]
    } else if (grepl("Reyfman2020|Angelidis2019", label)) {
      bulk_comparison = '1'
      bulk = bulk_de[[1]] %>% ungroup()
    } else {
      stop("not sure what to do with label: ", label)
    }
    
    # fix column names
    fix_colnames = function(df) {
      colnames(df) %<>%
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
      return(df)
    }
    sc %<>% fix_colnames()
    bulk %<>% fix_colnames()
    
    # call FPs
    ns_sc = sc %>%
      # replace Bonferroni with BH correction
      mutate(padj = p.adjust(p_val, 'BH')) %>%
      filter(padj > 0.1) %>%
      pull(gene)
    ns_bulk = filter(bulk, p_val_adj > 0.1) %>% pull(gene)
    
    # single-cell FPs
    fps = sc %>%
      arrange(p_val) %>%
      filter(gene %in% ns_bulk) %>%
      head(200) %>%
      mutate(rank = row_number(),
             sc_comparison = sc_comparison,
             bulk_comparison = bulk_comparison) %>%
      dplyr::select(sc_comparison, bulk_comparison, cell_type,
                    rank, gene, everything())
    if ("runtime" %in% colnames(fps)) {
      fps %<>% dplyr::select(-runtime, -mem_usage)
    }
    
    # single-cell FNs
    fn_genes = bulk %>%
      arrange(p_val) %>%
      filter(gene %in% ns_sc) %>%
      filter(!duplicated(gene)) %>%
      head(200) %>%
      pull(gene)
    fns = sc %>%
      filter(gene %in% fn_genes) %>%
      # order by bulk p-values
      right_join(data.frame(gene = fn_genes), by = 'gene') %>%
      mutate(rank = row_number(),
             sc_comparison = sc_comparison,
             bulk_comparison = bulk_comparison) %>%
      dplyr::select(sc_comparison, bulk_comparison, cell_type,
                    rank, gene, everything())
    if ("runtime" %in% colnames(fns)) {
      fns %<>% dplyr::select(-runtime, -mem_usage)
    }
    
    # merge in expression summary to both
    summary0 = filter(expr_summary,
                      cell_type == !!cell_type,
                      comparison == sc_comparison) %>%
      dplyr::rename(sc_comparison = comparison) %>%
      dplyr::select(-dataset)
    fps %<>% left_join(summary0, by = c('cell_type', 'sc_comparison', 'gene'))
    fns %<>% left_join(summary0, by = c('cell_type', 'sc_comparison', 'gene'))
    
    # append to results
    FP %<>% bind_rows(fps)
    FN %<>% bind_rows(fns)
  }
}

# construct output
output = list(FPs = FP, FNs = FN)

# create output directory, if it doesn't exist
output_dir = dirname(args$output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
# save results
saveRDS(output, args$output_file)
