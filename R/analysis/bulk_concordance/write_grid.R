# write concordance array for scRNA-seq/bulk comparisons
# FCC array
fcc_opts = list(
  method = 'fcc',
  cor_method = 'spearman'
)
fcc_array = do.call(expand.grid, c(fcc_opts, stringsAsFactors = F))
# AUCC array
aucc_opts = list(
  method = 'aucc',
  k = c(100, 200, 500, 1000)
)
aucc_array = do.call(expand.grid, c(aucc_opts, stringsAsFactors = F))
# create results template using all parameters of interest
template = bind_rows(
  fcc_array,
  aucc_array
)
