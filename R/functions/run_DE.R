## Run a DE analysis within each cell type using Seurat
## custom fork of Seurat is available at https://github.com/jordansquair/Seurat
seurat_DE = function(sc, de_test, params = list()) {
  library(Matrix)
  # get classes
  meta = sc@meta.data %>%
    # convert cell_type from factor to character, if applicable
    mutate(cell_type = as.character(cell_type))
  labels = unique(meta$label)
  if (is.factor(labels)) {
    ## https://github.com/satijalab/seurat/issues/741
    ## avg_logFC: log fold-chage of the average expression between the two
    ## groups. Positive values indicate that the gene is more highly expressed
    ## in the first group
    label1 = levels(labels)[2]
    label2 = levels(labels)[1]
  } else {
    label1 = labels[1]
    label2 = labels[2]
  }
  # get cell types
  cell_types = unique(meta$cell_type)
  # make sure idents are set right in the Seurat object
  Idents(sc) = meta$cell_type
  
  # check if integer or already normalized, normalize if needed
  mat = GetAssayData(sc, slot = 'counts')
  if ((sum(mat %% 1 == 0) == length(mat)) == T) {
    sc %<>% NormalizeData()
  } else {
    sc[['RNA']]@data = mat
  }
  
  # check if we are doing mixed model
  if (de_test %in% c("mixed_lm", "mixed_nbinom", "mixed_poisson")) {
    replicate.var = 'replicate'
  } else {
    replicate.var = NULL
  }
  
  # check for extra parameters for the mixed models
  if (!is.null(params$test)) {
    if (params$test == "LRT") {
      use.lrt = T
    } else {
      use.lrt = F
    }
  } else {
    use.lrt = F
  }
  
  if (!is.null(params$offset)) {
    if (params$offset == "YES") {
      use.offset = T
    } else {
      use.offset = F
    }
  } else {
    use.offset = F
  }
  
  # run Seurat DE
  DE = list()
  for (cell_type_idx in seq_along(cell_types)) {
    cell_type = cell_types[cell_type_idx]
    message("[", cell_type_idx, "/", length(cell_types),
            "] working on cell type: ", cell_type, " ...")
    
    # check to make sure there are enough cells
    n_cells = table(meta$cell_type, meta$label)
    if (min(n_cells[cell_type, ]) < 3) {
      message(" .. not enough cells, skipping ...")
      next
    } else {
      tryCatch({
        # subset to the right cell type
        Idents(sc) = sc$cell_type
        sub = sc %>% subset(idents = cell_type)
        # drop genes that are never expressed, as Seurat does
        keep = rowSums(sub) > 0
        sub = sub[keep, ]
        # run DE analysis
        markers = FindMarkers(sub, ident.1 = label1, ident.2 = label2,
                              assay = 'RNA', slot = 'data',
                              min.pct = -Inf, min.cells.feature = 0,
                              min.cells.group = 3, logfc.threshold = -Inf,
                              group.by = 'label', subset.ident = cell_type,
                              test.use = de_test, replicate.var = replicate.var,
                              use.offset = use.offset, use.lrt = use.lrt) %>%
          rownames_to_column('gene') %>%
          mutate(test = de_test)
        DE[[cell_type]] = markers
      }, error = function(e) message(e))
    }
  }
  return(DE)
}

## Run a DE analysis within each cell type using pseudobulk methods
pseudobulk_DE = function(sc,
                         de_test = c('pseudobulk_limma',
                                     'pseudobulk_DESeq2',
                                     'pseudobulk_edgeR'),
                         params = list()) {
  
  # metadata must contain replicate column
  meta = sc@meta.data %>%
    # convert cell_type from factor to character, if applicable
    mutate(cell_type = as.character(cell_type))
  if (!"replicate" %in% colnames(meta))
    stop("metadata does not contain replicate information")
  
  # keep only cell types with enough cells
  keep = meta %>%
    dplyr::count(cell_type, label) %>%
    group_by(cell_type) %>%
    filter(all(n >= 3)) %>%
    pull(cell_type) %>%
    unique()
  
  # process data into gene x replicate x cell_type matrices
  Idents(sc) = sc$cell_type
  pseudobulks = keep %>%
    map( ~ {
      print(.)
      sc_sub = subset(sc, idents = .)
      meta_sub = sc_sub@meta.data %>%
        mutate(label = as.character(label),
               replicate = as.character(replicate))
      # catch cell types without replicates or conditions
      if (n_distinct(meta_sub$label) < 2)
        return(NA)
      if (!is.null(params$replicate)) {
        # optionally, turn off replicate summarization
        if (params$replicate == 'cells') {
          meta_sub$replicate = colnames(sc_sub)
          # make sure there are no ":" in the rownames - causes a bug downstream
          if (any(grepl("\\:", meta_sub$replicate))) {
            meta_sub$replicate = gsub(".*\\:", "", meta_sub$replicate)
          }
        }
      } else {
        replicate_counts = distinct(meta_sub, label, replicate) %>%
          group_by(label) %>%
          summarise(replicates = n_distinct(replicate)) %>%
          pull(replicates)
        if (any(replicate_counts < 2))
          return(NA)
      }
      
      # process data into gene X replicate X cell_type matrices 
      # https://jef.works/blog/2020/04/06/quickly-creating-pseudobulks/
      mm = model.matrix(~ 0 + replicate:label, data = meta_sub)
      mat_mm = GetAssayData(sc_sub, slot = 'counts') %*% mm
      keep_genes = rowSums(mat_mm > 0) > 0
      mat_mm = mat_mm[keep_genes, ] %>% as.data.frame()
      mat_mm %<>% as.data.frame()
      colnames(mat_mm) = gsub("replicate|label", "", colnames(mat_mm))
      # drop empty columns
      keep_samples = colSums(mat_mm) > 0
      mat_mm %<>% extract(, keep_samples)
      return(mat_mm)
    }) %>%
    setNames(keep)
  # drop NAs
  pseudobulks %<>% extract(!is.na(.))
  
  # also filter out cell types with no retained genes
  min_dim = map(pseudobulks, as.data.frame) %>% map(nrow)
  pseudobulks %<>% extract(min_dim > 1)
  
  # also filter out types without replicates
  min_repl = map_int(pseudobulks, ~ {
    # make sure we have a data frame a not a vector
    tmp = as.data.frame(.)
    targets = data.frame(group_sample = colnames(tmp)) %>%
      mutate(group = gsub(".*\\:", "", group_sample))
    if (n_distinct(targets$group) == 1)
      return(as.integer(0))
    min(table(targets$group))
  })
  pseudobulks %<>% extract(min_repl >= 2)
  
  # run pseudobulk DE
  if (de_test == 'pseudobulk_limma') {
    library(limma)
    
    # run limma on each cell type
    DE = pseudobulks %>%
      map(function(x) {
        tryCatch({
          # create targets matrix
          targets = data.frame(group_sample = colnames(x)) %>%
            mutate(group = gsub(".*\\:", "", group_sample))
          ## optionally, carry over factor levels from entire dataset
          if (is.factor(meta$label)) {
            targets$group %<>% factor(levels = levels(meta$label))
          }
          if (n_distinct(targets$group) > 2)
            return(NULL)
          # create design
          design = model.matrix(~ group, data = targets)
          
          # limma-trend vs. limma-voom
          used_voom = FALSE
          if (params$mode == 'trend') {
            library(edgeR)
            dge = DGEList(as.matrix(x), group = targets$group)
            dge = calcNormFactors(dge)
            x = new("EList")
            x$E = edgeR::cpm(dge, log = TRUE, prior.count = 3)
          } else {
            ## limma-voom: default
            counts = all(as.matrix(x) %% 1 == 0)
            if (counts) {
              x = voom(as.matrix(x), design)
              used_voom = T
            }
          }
          # run lmFit
          trend_bool = params$mode == 'trend'
          fit = lmFit(x, design) %>%
            eBayes(trend = trend_bool, robust = trend_bool)
          # format the results
          res = fit %>%
            # extract all coefs except intercept
            topTable(number = Inf, coef = -1) %>%
            rownames_to_column('gene') %>%
            # flag voom usage
            mutate(used_voom = used_voom) %>%
            # flag test used
            mutate(test = 'pseudobulk_limma')
          return(res)
        }, error = function(e) {
          message(e)
          return(data.frame())
        })
      })
  } else if (de_test == 'pseudobulk_DESeq2') {
    library(DESeq2)
    
    # run DEseq2 for each cell_type
    DE = pseudobulks %>%
      map(function(x) {
        tryCatch({
          # create DESeq dataset
          targets = data.frame(group_sample = colnames(x)) %>%
            mutate(group = gsub(".*\\:", "", group_sample))
          ## optionally, carry over factor levels from entire dataset
          if (is.factor(meta$label)) {
            targets$group %<>% factor(levels = levels(meta$label))
          }
          dds = DESeqDataSetFromMatrix(countData = x,
                                       colData = targets,
                                       design = ~ group)
          
          # set defaults on params
          if (!'test' %in% names(params))
            params$test = 'LRT'
          if (!'fitType' %in% names(params))
            params$fitType = 'parametric'
          if (!'sfType' %in% names(params))
            params$sfType = 'poscounts'
          if (!'independentFiltering' %in% names(params))
            params$independentFiltering = FALSE
          if (!'betaPrior' %in% names(params))
            params$betaPrior = FALSE
          
          # run differential expression
          if (params$test == 'Wald') {
            # remove 'reduced' argument for Wald test
            dds = try(DESeq(dds,
                            test = params$test,
                            # reduced = ~ 1,
                            fitType = params$fitType,
                            sfType = params$sfType,
                            betaPrior = params$betaPrior))
            
          } else {
            dds = try(DESeq(dds,
                            test = params$test,
                            reduced = ~ 1,
                            fitType = params$fitType,
                            sfType = params$sfType,
                            betaPrior = params$betaPrior))
          }
          res = results(dds,
                        independentFiltering = params$independentFiltering)
          # write
          result = as.data.frame(res) %>%
            mutate(gene = rownames(x)) %>%
            # flag test used
            mutate(test = 'pseudobulk_DESeq2')
          return(result)
        }, error = function(e) {
          message(e)
          return(data.frame())
        })
      })
  } else if (de_test == 'pseudobulk_edgeR') {
    library(edgeR)
    # run edgeR for each cell_type
    DE = pseudobulks %>%
      map(function(x) {
        tryCatch({
          targets = data.frame(group_sample = colnames(x)) %>%
            mutate(group = gsub(".*\\:", "", group_sample))
          ## optionally, carry over factor levels from entire dataset
          if (is.factor(meta$label)) {
            targets$group %<>% factor(levels = levels(meta$label))
          }
          design = model.matrix(~ group, data = targets)
          y = DGEList(counts = x, group = targets$group)
          
          # set defaults on params
          if (!'test' %in% names(params))
            params$test = 'LRT'
          if (!'method' %in% names(params))
            params$method = 'TMM'
          if (!'trend.method' %in% names(params))
            params$trend.method = 'locfit'
          if (!'robust' %in% names(params))
            params$robust = FALSE
          if (!'tagwise' %in% names(params))
            params$tagwise = TRUE
          if (!'prior' %in% names(params))
            params$prior = NULL
          
          # catch some parameter type errors
          if (!is.null(params$prior)) {
            params$prior = as.numeric(params$prior)
          }
          if (params$tagwise == "FALSE") {
            params$tagwise = FALSE
          }
          
          y = calcNormFactors(y, method = params$method)
          if (params$robust) {
            y = estimateGLMRobustDisp(y, design,
                                      trend.method = params$trend.method)
          } else {
            y = estimateDisp(y, design, trend.method = params$trend.method,
                             tagwise = params$tagwise, prior.df = params$prior)
          }
          
          if (params$test == 'LRT') {
            fit = glmFit(y, design = design)
            test = glmLRT(fit)
          } else {
            # QLF: default
            fit = glmQLFit(y, design)
            test = glmQLFTest(fit, coef = -1)
          }
          res = topTags(test, n = Inf) %>%
            as.data.frame() %>%
            rownames_to_column('gene') %>%
            # flag test used
            mutate(test = 'pseudobulk_edgeR')
        }, error = function(e) {
          message(e)
          return(data.frame())
        })
      })
  }
  
  # standardize column names to match Seurat output
  DE %<>% map(~ {
    colnames(.) %<>%
      fct_recode('p_val' = 'p.value',  ## DESeq2
                 'p_val' = 'p.value',  ## t/wilcox
                 'p_val' = 'P.Value',  ## limma
                 'p_val' = 'PValue'  , ## edgeR
                 'p_val_adj' = 'padj', ## DESeq2/t/wilcox
                 'p_val_adj' = 'adj.P.Val',      ## limma
                 'p_val_adj' = 'FDR',            ## edgeER
                 'avg_logFC' = 'log2FoldChange', ## DESEeq2
                 'avg_logFC' = 'logFC' ## limma/edgeR
      ) %>%
      as.character()
    return(.)
  })
  
  return(DE)
}

## Run a DE analysis on a bulk dataset object including targets
bulk_DE = function(x, targets, de_test = c('bulk_limma',
                                           'bulk_DESeq2',
                                           'bulk_edgeR'),
                   used_voom = 'auto',
                   params = list()) {
  # split out params
  params_vec = strsplit(de_test, ',') %>%
    unlist() %>%
    extract(-1)
  for (param_idx in seq_along(params_vec)) {
    param_name = gsub("\\?.*$" ,"", params_vec[param_idx])
    param_value = gsub("^.*\\?" ,"", params_vec[param_idx])
    params[[param_name]] = param_value
  }
  de_test = gsub(",.*$", "", de_test)
  ## print them
  message("running ", de_test, " with params: ", names(params), params)
  
  if (de_test == 'bulk_limma') {
    library(limma)
    
    # run limma
    if (n_distinct(targets$label) > 2)
      return(NULL)
    # create design
    design = model.matrix(~ label, data = targets)
    
    # limma-trend vs. limma-voom
    used_voom = FALSE
    if (params$mode == 'trend') {
      library(edgeR)
      dge = DGEList(as.matrix(x), group = targets$label)
      dge = calcNormFactors(dge)
      x = new("EList")
      x$E = edgeR::cpm(dge, log = TRUE, prior.count = 3)
    } else {
      ## limma-voom: default
      counts = all(as.matrix(x) %% 1 == 0)
      if (counts) {
        x = voom(as.matrix(x), design)
        used_voom = T
      }
    }
    
    # run limma
    trend_bool = params$mode == 'trend'
    DE = lmFit(x, design) %>%
      eBayes(trend = trend_bool, robust = trend_bool) %>%
      # extract all coefs except intercept
      topTable(number = Inf, coef = -1) %>%
      rownames_to_column('gene') %>%
      # flag voom usage
      mutate(used_voom = used_voom) %>%
      # flag test used
      mutate(test = 'bulk_limma')
    return(DE)
  } else if (de_test == 'bulk_DESeq2') {
    library(DESeq2)
    
    # run DEseq2
    # create DESeq dataset
    dds = DESeqDataSetFromMatrix(countData = x,
                                 colData = targets,
                                 design = ~ label)
    
    # set defaults on params
    if (!'test' %in% names(params))
      params$test = 'LRT'
    if (!'fitType' %in% names(params))
      params$fitType = 'parametric'
    if (!'sfType' %in% names(params))
      params$sfType = 'poscounts'
    if (!'independentFiltering' %in% names(params))
      params$independentFiltering = FALSE
    if (!'betaPrior' %in% names(params))
      params$betaPrior = FALSE
    
    # run differential expression
    if (params$test == 'Wald') {
      # remove 'reduced' argument for Wald test
      dds = try(DESeq(dds,
                      test = params$test,
                      # reduced = ~ 1,
                      fitType = params$fitType,
                      sfType = params$sfType,
                      betaPrior = params$betaPrior))
      
    } else {
      dds = try(DESeq(dds,
                      test = params$test,
                      reduced = ~ 1,
                      fitType = params$fitType,
                      sfType = params$sfType,
                      betaPrior = params$betaPrior))
    }
    res = results(dds,
                  independentFiltering = params$independentFiltering)
    
    # write
    DE = as.data.frame(res) %>%
      mutate(gene = rownames(x)) %>%
      # flag test used
      mutate(test = 'bulk_DESeq2')
    return(DE)
  } else if (de_test == 'bulk_edgeR') {
    library(edgeR)
    
    # run edgeR
    design = model.matrix(~ label, data = targets)
    y = DGEList(counts = x, group = targets$label)
    
    # set defaults on params
    if (!'test' %in% names(params))
      params$test = 'LRT'
    if (!'method' %in% names(params))
      params$method = 'TMM'
    if (!'trend.method' %in% names(params))
      params$trend.method = 'locfit'
    if (!'robust' %in% names(params))
      params$robust = FALSE
    
    y = calcNormFactors(y, method = params$method)
    if (params$robust) {
      y = estimateGLMRobustDisp(y, design,
                                trend.method = params$trend.method)
    } else {
      y = estimateDisp(y, design, trend.method = params$trend.method)
    }
    
    if (params$test == 'LRT') {
      fit = glmFit(y, design = design)
      test = glmLRT(fit)
    } else {
      # QLF: default
      fit = glmQLFit(y, design)
      test = glmQLFTest(fit, coef = -1)
    }
    
    # write
    DE = topTags(test, n = Inf) %>%
      as.data.frame() %>%
      rownames_to_column('gene') %>%
      # flag test used
      mutate(test = 'bulk_edgeR')
    return(DE)
  }
  
  # standardize column names to match Seurat output
  colnames(DE) %<>%
    fct_recode('p_val' = 'p.value',  ## DESeq2
               'p_val' = 'p.value',  ## t/wilcox
               'p_val' = 'P.Value',  ## limma
               'p_val' = 'PValue'  , ## edgeR
               'p_val_adj' = 'padj', ## DESeq2/t/wilcox
               'p_val_adj' = 'adj.P.Val',      ## limma
               'p_val_adj' = 'FDR',            ## edgeER
               'avg_logFC' = 'log2FoldChange', ## DESEeq2
               'avg_logFC' = 'logFC' ## limma/edgeR
    ) %>%
    as.character()
  return(DE)
}

## Run a scRNAseq differential expression analysis, given a Seurat object.
run_DE = function(sc, de_test, params = list(), ...) {
  library(peakRAM)
  # check there are exactly two labels
  meta = sc@meta.data
  n_labels = n_distinct(meta$label)
  if (n_labels != 2) {
    labels = paste0(unique(meta$label), collapse = ' | ')
    stop("can't run pairwise DE on ", n_labels, " labels (", labels, ")")
  }
  
  # split out params
  params_vec = strsplit(de_test, ',') %>%
    unlist() %>%
    extract(-1)
  for (param_idx in seq_along(params_vec)) {
    param_name = gsub("\\?.*$" ,"", params_vec[param_idx])
    param_value = gsub("^.*\\?" ,"", params_vec[param_idx])
    params[[param_name]] = param_value
  }
  
  de_test = gsub(",.*$", "", de_test)
  ## print them
  message("running ", de_test, " with params: ", names(params), params)
  
  # run DE tests
  if (grepl("pseudobulk", de_test)) {
    time = system.time({
      mem = peakRAM({
        DE = pseudobulk_DE(sc, de_test, params) %>%
          bind_rows(.id = 'cell_type')
      })
    })
    DE %<>% mutate(runtime = time[3], mem_usage = mem[1, 4])
  } else {
    DE = data.frame()
    tryCatch({
      time = system.time({
        mem = peakRAM({
          DE = seurat_DE(sc, de_test, params) %>%
            bind_rows(.id = 'cell_type')
        })
      })
      DE %<>% mutate(runtime = time[3], mem_usage = mem[1, 4])
    }, error = function(e) message(e))
  }
  return(DE)
}
