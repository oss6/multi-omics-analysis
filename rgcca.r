library(RGCCA)

source('./fileio.R')
source('./process_peaks.R')

run_rgcca <- function(genes, response, use_specific_peaks, results_subdir = '') {
  rna.ds <- load.dataset(
    meta.file = './data/sample_sheet.csv', meta.sep = ',',
    data.file = './data/rna_norm_counts.csv', data.sep = ','
  )
  
  # rna.ds$meta.data <- rna.ds$meta.data[!grepl('C', rownames(rna.ds$meta.data)),]
  # rna.ds$data.matrix <- rna.ds$data.matrix[!grepl('C', rownames(rna.ds$data.matrix)),]
  
  ppolar.ds <- load.dataset(
    meta.file = './data/sample_sheet.csv', meta.sep = ',',
    data.file = './data/polar_pos_pqn_imputed_glog.csv', data.sep = ','
  )
  
  # ppolar.ds$meta.data <- ppolar.ds$meta.data[!grepl('B', rownames(ppolar.ds$meta.data)),]
  # ppolar.ds$data.matrix <- ppolar.ds$data.matrix[!grepl('B', rownames(ppolar.ds$data.matrix)),]
  
  site <- data.matrix(data.frame(rna.ds$meta.data$Site))
  rownames(site) <- row.names(rna.ds$meta.data)
  colnames(site) <- c('site')
  
  ref <- data.matrix(data.frame(rna.ds$meta.data$REF))
  rownames(ref) <- row.names(rna.ds$meta.data)
  colnames(ref) <- c('ref')
  
  # wc1 <- rep(t(water_chemicals[1,]), each = 6)
  # wc1 <- data.matrix(c(wc1, wc1))
  # rownames(wc1) <- row.names(rna.ds$meta.data)
  # colnames(wc1) <- c('wc1')
  
  responses <- list(
    site = site,
    ref = ref,
    wc1 = wc1
  )
  response_map <- list(
    site = 'Site',
    ref = 'REF',
    wc1 = 'wc1'
  )
  
  shared.ids <- intersect(row.names(rna.ds$data.matrix), row.names(ppolar.ds$data.matrix))
  
  # Genes blocks
  gene_blocks <- list()
  
  if (is.list(genes)) {
    gene_blocks <- lapply(genes, function (x) { rna.ds$data.matrix[shared.ids, colnames(rna.ds$data.matrix) %in% x] })
  } else {
    gene_blocks <- list(rna = rna.ds$data.matrix[shared.ids, colnames(rna.ds$data.matrix) %in% genes])
  }
  
  # Peaks blocks
  peak_blocks <- list()
  
  if (use_specific_peaks) {
    peak_blocks <- lapply(get_partitioned_peaks(), function (x) {
      ppolar.ds$data.matrix[shared.ids, colnames(ppolar.ds$data.matrix) %in% x]
    })
  } else {
    peak_blocks[['ppolar']] = ppolar.ds$data.matrix[shared.ids,]
  }
  
  # Response block
  response_block <- list()
  response_block[[response]] <- data.matrix(responses[[response]][shared.ids, ])
  response_list <- rna.ds$meta.data[shared.ids,][[response_map[[response]]]]
  
  blocks <- c(
    gene_blocks,
    peak_blocks,
    response_block
  )
  
  # cca_perm_res <- rgcca_permutation(
  #   blocks = A,
  #   connection = C,
  #   superblock = FALSE,
  #   response = length(A),
  #   par_type = 'sparsity',
  #   scheme = 'factorial',
  #   method = 'sgcca',
  #   n_cores = 6
  # )
  
  ncomp <- c(rep(4, length(blocks) - 1), 1)
  
  cca_cv <- rgcca_cv(
    blocks = blocks,
    response = length(blocks),
    ncomp = ncomp,
    k = 20,
    method = 'sgcca',
    par_type = 'sparsity',
    scheme = 'factorial',
    n_cores = 6
  )
  
  cca.res <- rgcca(
    blocks = blocks,
    method = 'sgcca',
    sparsity = cca_cv$bestpenalties,
    scheme = 'factorial',
    ncomp = ncomp,
    response = length(blocks),
    superblock = FALSE,
    verbose = FALSE)
  
  # Save results
  results_dir <- paste('results', 'rgcca', results_subdir, sep = '/')
  if (!dir.exists(results_dir)) {
    dir.create(results_dir)
  }
  
  prefix = paste(response, ifelse(use_specific_peaks, 'specpeaks', 'allpeaks'), sep = '_')
  
  filename = paste(paste(prefix, 'ave', sep = '_'), 'png', sep = '.')
  png(file=paste(results_dir, filename, sep = '/'))
  plot(cca.res, type = 'ave')
  dev.off()
  
  if (use_specific_peaks) {
    for (block in seq(2, length(blocks) - 1)) {
      filename = paste(paste(prefix, 'components', 'rna', 'vs', names(blocks)[block], sep = '_'), 'png', sep = '.')
      png(file=paste(results_dir, filename, sep = '/'))
      plot(cca.res, type = "samples", comp = 1, block = c(1, block), response = response_list)
      dev.off()
    }
  } else {
    for (comp in c(1, 2, 3)) {
      filename = paste(paste(prefix, 'components', comp, sep = '_'), 'png', sep = '.')
      png(file=paste(results_dir, filename, sep = '/'))
      plot(cca.res, type = "samples", comp = comp, block = c(1, 2), response = response_list)
      dev.off()
    }
  }
  
  # stability_selection <- rgcca_stability(degs_with_metabolites_cca)
  boot_res <- rgcca_bootstrap(cca.res, n_cores = 6, verbose = FALSE)
  
  return(list(
    sgcca_res=cca.res,
    boot_res=boot_res
  ))
}

