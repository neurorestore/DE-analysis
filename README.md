# DE-analysis

This repository contains R source code used to conduct the analysis in our manuscript, "Confronting false discoveries in single-cell differential expression."

A brief overview of the main computational analyses that were conducted, and the location of the corresponding source code, is given below.

- First, differentially expressed genes were identified in single-cell and matching bulk datasets, respectively, using code in the directories `R/analysis/run_de` and `R/analysis/run_bulk_de`. Code in `R/analysis/run_spike_in_de` was used to analyze a lone single-cell dataset in which the ERCC mixture of synthetic mRNAs was spiked in alongside each individual cell.
    - A list of all single-cell and bulk datasets analyzed in this study is provided in `R/functions/datasets.R`. Datasets containing multiple comparisons of two experimental groups were split into each possible comparison using the functions in `R/functions/get_comparisons.R` and `R/functions/get_bulk_comparisons.R`. 
    - Code used to run differential expression analyses is provided in `R/functions/run_DE.R`
- The concordance between the single-cell and bulk DE results was then quantified using code in the `R/analysis/bulk_concordance` directory.
    - Code used to calculate the AUCC and fold-change correlation is provided in `R/functions/calculate_overlap.R`
- Gene set enrichment analysis of differential expression results was performed for both the single-cell and bulk datasets using code in the `R/analysis/run_GSEA` directory.
- False-positive and false-negative DE calls were obtained for the single-cell data, using the bulk data as a reference, using code in the `R/analysis/extract_FPs` directory.
- A number of summary statistics were obtained for each dataset (e.g., number of replicates, number of cell types) or each gene within each dataset (e.g., mean expression, delta-variance), using code in the directories `R/analysis/confounds`, `R/analysis/delta_variance`, and `R/analysis/expr_summary`. 
- The relationships between mean expression, the variance of gene expression, and the delta-variance in 'pseudo-replicates' were interrogated using code in `R/analysis/mean_variance`. 
- The effect of between-replicate variance was interrogated with simulation studies using the code in `R/analysis/simulations`. The code in this directory was used to generate synthetic gene expression data, perform DE analysis, and analyse the properties of DE genes. 
- DE analysis was performed between random groups of control samples using code in the `R/analysis/control_only` directory. This code was also used to analyze a spatial transcriptomics dataset.
- The performance of generalized linear mixed models was assessed in downsampled datasets using code in the `R/analysis/downsample_cells` directory.
- Finally, the computational resources (wall time, peak RAM usage) used by each method were extracted using code in the `R/analysis/time_RAM` directory.
