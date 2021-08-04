get_bulk_comparisons = function(dataset, expr, meta) {
  # set up container
  comparisons = list()
  # handle each dataset appropriately
  if (dataset == 'Angelidis2019') {
    meta %<>% mutate(label = factor(label, levels = c("3m", "24m")))
    results[[1]] = list(expr = expr, meta = meta)
  } else if (dataset == 'Hagai2018') {
    ## Hagai2018: two different binary comparisons
    for (comparison in c('LPS4', 'PIC4')) {
      message('  processing comparison: ', comparison, ' ...')
      meta0 = meta %>%
        rownames_to_column(var = 'sample') %>%
        filter(label %in% c('UNST', comparison)) %>%
        mutate(label = factor(label, levels = c("UNST", comparison)))
      expr0 = expr[, meta0$sample]
      results[[comparison]] = list(expr = expr0, meta = meta0)
    }
  } else if (dataset == 'CanoGamez2020') {
    ## CanoGamez2020: Two different cell types, 7 different cytokine conditions
    cytokines = c("IFNB", "Th17", "Resting", "Th2", "Th0", "iTreg", "Th1")
    cell_types = c('Naive', 'Memory')
    stimulation_times = c("16h", "5d")
    grid = tidyr::crossing(cytokine1 = cytokines,
                           cytokine2 = cytokines,
                           cell_type = cell_types,
                           stimulation_time = stimulation_times) %>%
      filter(cytokine1 != cytokine2) %>%
      filter(cytokine1 == 'Resting')
    for (grid_idx in seq_len(nrow(grid))) {
      cytokine1 = grid$cytokine1[grid_idx]
      cytokine2 = grid$cytokine2[grid_idx]
      cell_type = grid$cell_type[grid_idx]
      stimulation_time = grid$stimulation_time[grid_idx]
      key = paste0(cytokine1, '|', cytokine2, "|", cell_type, "|",
                   stimulation_time)
      message('  processing comparison: ', key, ' ...')
      
      meta0 = meta %>%
        mutate(idx = row_number()) %>%
        filter(stimulation_time == !!stimulation_time) %>%
        filter(grepl(!!cell_type, cell_type)) %>%
        filter(cytokine_condition %in% c(cytokine1, cytokine2)) %>%
        mutate(cytokine_condition = factor(cytokine_condition,
                                           levels = c(cytokine1, cytokine2)),
               label = cytokine_condition)
      expr0 = expr %>% extract(, meta0$idx)
      results[[key]] = list(expr = expr0, meta = meta0)
    }
  } else if (dataset == 'CanoGamez2020:proteomics') {
    ## CanoGamez2020: Two different cell types, 7 different cytokine conditions
    cytokines = c("IFNB", "Th17", "Resting", "Th2", "Th0", "iTreg", "Th1")
    cell_types = c('Naive', 'Memory')
    grid = tidyr::crossing(cytokine1 = cytokines,
                           cytokine2 = cytokines,
                           cell_type = cell_types) %>%
      filter(cytokine1 != cytokine2) %>%
      filter(cytokine1 == 'Resting')
    for (grid_idx in seq_len(nrow(grid))) {
      cell_type = grid$cell_type[grid_idx]
      cytokine1 = grid$cytokine1[grid_idx]
      cytokine2 = grid$cytokine2[grid_idx]
      key = paste0(cytokine1, '|', cytokine2, "|", cell_type)
      message('  processing comparison: ', key, ' ...')
      
      meta0 = meta %>%
        mutate(idx = row_number()) %>%
        filter(cell_type == !!cell_type) %>%
        filter(cytokine_condition %in% c(cytokine1, cytokine2)) %>%
        mutate(cytokine_condition = factor(cytokine_condition,
                                           levels = c(cytokine1, cytokine2)),
               label = cytokine_condition)
      expr0 = expr %>% extract(, meta0$idx)
      results[[key]] = list(expr = expr0, meta = meta0)
    }
  } else {
    stop("invalid dataset: ", dataset, " ...")
  }
  
  # drop all unused factor levels
  for (comparison_idx in seq_along(results)) {
    results[[comparison_idx]]$meta %<>% droplevels()
  }
  
  return(results)
}
