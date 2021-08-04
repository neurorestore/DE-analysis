## Get subset expression matrices containing all comparisons for a given
## dataset.
get_comparisons = function(dataset, expr, meta) {
  # set up container
  results = list()
  # handle each dataset appropriately
  if (dataset %in% c('Arneson2018',
                     'Avey2018',
                     'Brenner2020',
                     'Cheng2019',
                     'Co2020',
                     'Crowell2019',
                     'Der2019_kidney',
                     'Der2019_skin',
                     'Grubman2019',
                     'Hashimoto2019',
                     'Hu2017',
                     'Jakel2019',
                     'Mathys2019',
                     'Nagy2020',
                     'OrdovasMontanes2018',
                     'Rault2020',
                     'Rossi2019',
                     'Sathyamurthy2018',
                     'Schafflick2020_CSF',
                     'Schafflick2020_PBMCs',
                     'Skinnider2020',
                     'Wang2020'
  )) {
    results[[1]] = list(expr = expr, meta = meta)
  } else if (dataset %in% c('Goldfarbmuren2020',
                            'Schirmer2019',
                            'Ximerakis2019')) {
    ## two cell type levels
    results[['cell_type']] = list(expr = expr, meta = meta)
    meta0 = meta %>%
      dplyr::select(-cell_type) %>%
      dplyr::rename(cell_type = global_cell_type)
    results[['global_cell_type']] = list(expr = expr, meta = meta0)
  } else if (dataset == 'Bhattacherjee2019') {
    ## Bhattacherjee2019: two possible levels of cell types, and
    ## three different timepoints
    timepoints = c('Maintenance', '48h', '15d')
    cell_types = c('cell_type', 'global_cell_type')
    grid = tidyr::crossing(timepoint = timepoints, cell_type = cell_types)
    for (grid_idx in seq_len(nrow(grid))) {
      timepoint = grid$timepoint[grid_idx]
      cell_type = grid$cell_type[grid_idx]
      key = paste0(timepoint, '|', cell_type)
      message('  processing comparison: ', key, ' ...')

      meta0 = meta %>%
        mutate(idx = row_number()) %>%
        filter(grepl(timepoint, label))
      expr0 = expr %>% extract(, meta0$idx)
      if (cell_type == "global_cell_type") {
        meta0 %<>%
          dplyr::select(-cell_type) %>%
          dplyr::rename(cell_type = global_cell_type)
      }

      results[[key]] = list(expr = expr0, meta = meta0)
    }
  } else if (dataset == 'Huang2020') {
    ## Huang2020: two possible levels of cell types, and
    ## three different comparisons
    conditions = c('CD', 'colitis', 'UC')
    cell_types = c('cell_type', 'global_cell_type')
    grid = tidyr::crossing(condition = conditions, cell_type = cell_types)
    for (grid_idx in seq_len(nrow(grid))) {
      condition = grid$condition[grid_idx]
      cell_type = grid$cell_type[grid_idx]
      key = paste0(condition, '|', cell_type)
      message('  processing comparison: ', key, ' ...')

      meta0 = meta %>%
        mutate(idx = row_number()) %>%
        filter(label %in% c(condition, 'control'))
      expr0 = expr %>% extract(, meta0$idx)
      if (cell_type == "global_cell_type") {
        meta0 %<>%
          dplyr::select(-cell_type) %>%
          dplyr::rename(cell_type = global_cell_type)
      }

      results[[key]] = list(expr = expr0, meta = meta0)
    }
  } else if (dataset == 'Reyes2020') {
    ## Reyes2020: two possible levels of cell types, and
    ## three different comparisons
    cohorts = c('ICU-SEP vs. ICU-NoSEP',
                'Sepsis vs. control',
                'Sepsis vs. Leuk-UTI')
    cell_types = c('cell_type', 'global_cell_type')
    grid = tidyr::crossing(cohort = cohorts, cell_type = cell_types)
    for (grid_idx in seq_len(nrow(grid))) {
      cohort = grid$cohort[grid_idx]
      cell_type = grid$cell_type[grid_idx]
      key = paste0(cohort, '|', cell_type)
      message('  processing comparison: ', key, ' ...')

      if (cohort == 'ICU-SEP vs. ICU-NoSEP') {
        meta0 = meta %>%
          mutate(idx = row_number()) %>%
          filter(label %in% c("ICU-SEP", "ICU-NoSEP"))
      } else if (cohort == 'Sepsis vs. control') {
        meta0 = meta %>%
          mutate(idx = row_number()) %>%
          filter(label %in% c("Int-URO", "URO", "Bac-SEP", "ICU-SEP",
                              "Control")) %>%
          mutate(label = ifelse(label == 'Control', label, 'Sepsis'))
      } else if (cohort == 'Sepsis vs. Leuk-UTI') {
        meta0 = meta %>%
          mutate(idx = row_number()) %>%
          filter(label %in% c("Int-URO", "URO", "Bac-SEP", "ICU-SEP",
                              "Leuk-UTI")) %>%
          mutate(label = ifelse(label == 'Leuk-UTI', label, 'Sepsis'))
      }
      expr0 = expr %>% extract(, meta0$idx)
      if (cell_type == "global_cell_type") {
        meta0 %<>%
          dplyr::select(-cell_type) %>%
          dplyr::rename(cell_type = global_cell_type)
      }
      results[[key]] = list(expr = expr0, meta = meta0)
    }
  } else if (dataset == 'Wu2017') {
    ## Wu2017: two different binary comparisons
    for (comparison in c('stress', 'seizure')) {
      message('  processing comparison: ', comparison, ' ...')
      meta0 = meta %>%
        mutate(idx = row_number()) %>%
        filter(label %in% c('control', comparison)) %>%
        droplevels()
      expr0 = expr[, meta0$idx]
      results[[comparison]] = list(expr = expr0, meta = meta0)
    }
  } else if (dataset == 'Wagner2018') {
    ## Wagner2018: tyrosinase vs. chordin
    meta0 = meta %>%
      mutate(idx = row_number()) %>%
      filter(label != 'WT')
    expr0 = expr[, meta0$idx]
    results[[1]] = list(expr = expr0, meta = meta0)
  } else if (dataset == 'Gunner2019') {
    ## Gunner2019: control vs. lesion, in entire dataset or by genotype
    comparisons = c('entire_dataset', 'homozygous', 'heterozygous')
    for (comparison in comparisons) {
      message('  processing comparison: ', comparison, ' ...')
      if (comparison == 'entire_dataset') {
        meta0 = meta
        expr0 = expr
      } else {
        meta0 = meta %>%
          mutate(idx = row_number()) %>%
          filter(genotype == Hmisc::capitalize(comparison))
        expr0 = expr[, meta0$idx]
      }
      results[[comparison]] = list(expr = expr0, meta = meta0)
    }
  } else if (dataset == 'Haber2017_droplet') {
    ## Haber2017: three binary comparisons vs. control
    comparisons = c('Hpoly.Day3', 'Hpoly.Day10', 'Salmonella')
    for (comparison in comparisons) {
      message('  processing comparison: ', comparison, ' ...')
      meta0 = meta %>%
        mutate(idx = row_number()) %>%
        filter(label %in% c('Control', comparison))
      expr0 = expr[, meta0$idx]
      results[[comparison]] = list(expr = expr0, meta = meta0)
    }
  } else if (dataset == 'Aztekin2019') {
    ## Aztekin2019: amputation response
    grps = c("ST40_1", "ST40_0")
    meta0 = meta %>%
      mutate(idx = row_number()) %>%
      filter(label %in% grps)
    expr0 = expr %>%
      extract(, meta0$idx)
    results[[1]] = list(expr = expr0, meta = meta0)
  } else if (dataset == 'Kim2019') {
    ## Kim2019: aggression
    meta0 = meta %>%
      mutate(idx = row_number()) %>%
      filter(label %in% c('Control', 'Aggression'))
    expr0 = expr[, meta0$idx]
    results[[1]] = list(expr = expr0, meta = meta0)
  } else if (dataset == 'Wirka2019') {
    ## Wirka2019: 8w/0w in WT
    # subset metadata
    meta0 = meta %>%
      mutate(idx = row_number())
    ## filter by genotype
    meta0 %<>% filter(phenotype == 'wt')
    ## filter by timepoint
    meta0 %<>% filter(label != '16wk')
    # subset expression
    expr0 = expr[, meta0$idx]
    results[[1]] = list(expr = expr0, meta = meta0)
  } else if (dataset == 'Jaitin2018_HFD') {
    ## Jaitin2019 (dataset 1): HFD vs. NC, 6w
    meta0 = meta %>%
      mutate(idx = row_number()) %>%
      filter(timepoint == 6)
    expr0 = expr[, meta0$idx]
    results[[1]] = list(expr = expr0, meta = meta0)
  } else if (dataset == 'CanoGamez2020') {
    ## CanoGamez2020: compare all cytokines to unstimulated
    comparisons = unique(meta$label) %>%
      setdiff('UNS')
    for (comparison in comparisons) {
      message('  processing comparison: ', comparison, ' ...')
      meta0 = meta %>%
        mutate(idx = row_number()) %>%
        filter(label %in% c('UNS', comparison)) %>%
        mutate(label = factor(label, levels = c('UNS', comparison)))
      expr0 = expr[, meta0$idx]
      results[[comparison]] = list(expr = expr0, meta = meta0)
    }
  } else if (dataset == 'Davie2018') {
    ## Davie2018: all combinations of age/sex/genotype
    ages = unique(meta$label)

    ### age
    comparisons = tidyr::crossing(age1 = ages, age2 = ages) %>%
      filter(age1 < age2)
    for (grid_idx in seq_len(nrow(comparisons))) {
      age1 = comparisons$age1[grid_idx]
      age2 = comparisons$age2[grid_idx]
      comparison = paste0(age1, '|', age2)
      message('  processing comparison: ', comparison, ' ...')
      meta0 = meta %>%
        mutate(idx = row_number()) %>%
        filter(label %in% c(age1, age2))
      expr0 = expr[, meta0$idx]
      results[[comparison]] = list(expr = expr0, meta = meta0)
    }

    ### sex
    message('  processing comparison: sex ...')
    meta0 = meta %>%
      dplyr::select(-label) %>%
      dplyr::rename(label = gender)
    results[['sex']] = list(expr = expr, meta = meta0)

    ### genotype
    message('  processing comparison: genotype ...')
    meta0 = meta %>%
      dplyr::select(-label) %>%
      dplyr::rename(label = genotype)
    results[['genotype']] = list(expr = expr, meta = meta0)
  } else if (dataset == 'Wagner2018') {
    ## Wagner2018: tyrosinase vs. chordin
    meta0 = meta %>%
      mutate(idx = row_number()) %>%
      filter(label != 'WT')
    expr0 = expr[, meta0$idx]
    results[[1]] = list(expr = expr0, meta = meta0)
  } else if (dataset == 'Hrvatin2018') {
    ## Hrvatin2018: focus on 0h vs. 4h
    meta0 = meta %>%
      mutate(idx = row_number()) %>%
      filter(label != '1h')
    expr0 = expr[, meta0$idx]
    results[[1]] = list(expr = expr0, meta = meta0)
  } else if (dataset == 'Madissoon2020') {
    ## Madissoon2020: compare all timepoints to 0 h
    comparisons = unique(meta$label) %>%
      setdiff('0h')
    for (comparison in comparisons) {
      message('  processing comparison: ', comparison, ' ...')
      meta0 = meta %>%
        mutate(idx = row_number()) %>%
        filter(label %in% c('0h', comparison)) %>%
        mutate(label = factor(label, levels = c('0h', comparison)))
      expr0 = expr[, meta0$idx]
      results[[comparison]] = list(expr = expr0, meta = meta0)
    }
  } else if (dataset == 'Tran2019') {
    ## Tran2019: compare all timepoints to control
    comparisons = unique(meta$label) %>% setdiff('Ctrl')
    for (comparison in comparisons) {
      message('  processing comparison: ', comparison, ' ...')
      meta0 = meta %>%
        mutate(idx = row_number()) %>%
        filter(label %in% c('Ctrl', comparison))
      expr0 = expr[, meta0$idx]
      results[[comparison]] = list(expr = expr0, meta = meta0)
    }
  } else if (dataset == 'Cuomo2020') {
    ## Cuomo2020: compare all timepoints to day 0
    comparisons = unique(meta$label) %>% setdiff('day0')
    for (comparison in comparisons) {
      message('  processing comparison: ', comparison, ' ...')
      meta0 = meta %>%
        mutate(idx = row_number()) %>%
        filter(label %in% c('day0', comparison)) %>%
        mutate(label = factor(label, levels = c('day0', comparison)))
      expr0 = expr[, meta0$idx]
      results[[comparison]] = list(expr = expr0, meta = meta0)
    }
  } else if (dataset %in% c(
    "Hagai2018_mouse",
    "Hagai2018_rat",
    "Hagai2018_pig",
    "Hagai2018_rabbit"
  )) {
    ## Hagai2018: compare all timepoints to unstimulated
    comparisons = unique(meta$label) %>% setdiff('unst')
    for (comparison in comparisons) {
      message('  processing comparison: ', comparison, ' ...')
      meta0 = meta %>%
        mutate(idx = row_number()) %>%
        filter(label %in% c('unst', comparison)) %>%
        mutate(label = factor(label, levels = c('unst', comparison)))
      expr0 = expr[, meta0$idx]
      results[[comparison]] = list(expr = expr0, meta = meta0)
    }
  } else if (dataset == 'Wilk2020') {
    ## Wilk2020: Compare each vent status to control, all COVID to control
    # and vent to no vent covid
    ## Comparison 1: all COVID to healthy control
    ## Comparison 2: ventilated COVID to healthy control
    ## Comparison 3: non-ventilated COVID to healthy control
    ## Comparison 4: non-ventilated COVID to ventilated covid
    comparisons = c("Healthy_COVID", "Healthy_vCOVID", "Healthy_nvCOVID",
                    "nvCOVID_vCOVID")
    for (comparison in comparisons) {
      message('  processing comparison: ', comparison, ' ...')
      if (comparison == "Healthy_COVID") {
        meta0 = meta %>%
          mutate(idx = row_number()) %>%
          mutate(label = factor(label, levels = c('Healthy', 'COVID')))
        expr0 = expr[, meta0$idx]
        results[[comparison]] = list(expr = expr0, meta = meta0)
      }
      if (comparison == "Healthy_vCOVID") {
        meta0 = meta %>%
          mutate(idx = row_number()) %>%
          filter(Ventilated %in% c("Healthy", "Vent")) %>%
          mutate(label = factor(label, levels = c('Healthy', 'COVID')))
        expr0 = expr[, meta0$idx]
        results[[comparison]] = list(expr = expr0, meta = meta0)
      }
      if (comparison == "Healthy_nvCOVID") {
        meta0 = meta %>%
          mutate(idx = row_number()) %>%
          filter(Ventilated %in% c("Healthy", "NonVent")) %>%
          mutate(label = factor(label, levels = c('Healthy', 'COVID')))
        expr0 = expr[, meta0$idx]
        results[[comparison]] = list(expr = expr0, meta = meta0)
      }
      if (comparison == "nvCOVID_vCOVID") {
        meta0 = meta %>%
          mutate(idx = row_number()) %>%
          filter(label == 'COVID') %>%
          mutate(label = Ventilated) %>%
          mutate(label = factor(label, levels = c('NonVent', 'Vent')))
        expr0 = expr[, meta0$idx]
        results[[comparison]] = list(expr = expr0, meta = meta0)
      }
    }
  } else if (dataset %in% c('Kotliarov2020')) {
    ## explicitly specify factor levels
    meta %<>% mutate(label = factor(label, levels = c('low', 'high')))
    results[[1]] = list(expr = expr, meta = meta)
  } else if (dataset == 'Kang2018') {
    ## explicitly specify factor levels
    meta %<>% mutate(label = factor(label, levels = c('ctrl', 'stim')))
    results[[1]] = list(expr = expr, meta = meta)
  } else if (dataset == 'Angelidis2019') {
    ## explicitly specify factor levels
    meta %<>% mutate(label = factor(label, levels = c('3m', '24m')))
    results[[1]] = list(expr = expr, meta = meta)
  } else if (dataset == 'Reyfman2020') {
    ## explicitly specify factor levels
    meta %<>% mutate(label = factor(label, levels = c('Control',
                                                      'Pulmonary fibrosis')))
    results[[1]] = list(expr = expr, meta = meta)
  } else if (dataset == 'Denisenko2020') {
    #' warm vs. cold (fresh): for comparison to bulk
    #' methanol vs. fresh, cryopreserved vs. fresh in warm/cold
    comparisons = c("warm_vs_cold",
                    "methanol_warm",
                    "methanol_cold",
                    "cryopreserved_warm",
                    "cryopreserved_cold")
    for (comparison in comparisons) {
      meta0 = meta %>%
        mutate(idx = row_number()) %>%
        # drop single-nucleus, v2/v3 comparison
        filter(!grepl("^SN|^SC", label))
      if (grepl("methanol", comparison)) {
        temperature = gsub("^.*_", "", comparison)
        meta0 %<>%
          filter(is.na(label2) | label2 == 'MeOH') %>%
          replace_na(list(label2 = 'fresh')) %>%
          filter(label == temperature) %>%
          mutate(label = factor(label2, levels = c('fresh', 'MeOH')))
      } else if (grepl("cryo", comparison)) {
        temperature = gsub("^.*_", "", comparison)
        meta0 %<>%
          filter(is.na(label2) | label2 == 'Cryo') %>%
          replace_na(list(label2 = 'fresh')) %>%
          filter(label == temperature) %>%
          mutate(label = factor(label2, levels = c('fresh', 'Cryo')))
      } else {
        # warm v. cold
        meta0 %<>%
          filter(is.na(label2)) %>%
          # set factor levels
          mutate(label = factor(label, levels = c('cold', 'warm')))
      }
      expr0 = expr[, meta0$idx]
      results[[comparison]] = list(expr = expr0, meta = meta0)
    }
    ##
  } else if (dataset == "Hagai2018_plate") {
    ## explicitly specify factor levels
    meta %<>% mutate(label = factor(label, levels = c(2, 6)))
    results[[1]] = list(expr = expr, meta = meta)
  } else if (dataset == "Maniatis2019_mouse") { 
    meta0 = meta %>%
      mutate(idx = row_number()) %>%
      # rename replicate, label, and region for compatibility
      dplyr::rename(replicate = isolate, label = breed, cell_type = region)
    results[[1]] = list(expr = expr, meta = meta0)
  } else {
    stop("invalid dataset: ", dataset, " ...")
  }

  # drop all unused factor levels
  for (comparison_idx in seq_along(results)) {
    results[[comparison_idx]]$meta %<>% droplevels()
  }

  return(results)
}
