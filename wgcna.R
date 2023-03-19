library(WGCNA)
library(dplyr)

source('./fileio.R')
source('./id_mapping.R')
source('./genes_annotation.R')

run_wgcna <- function () {
  # allow multi-threading within WGCNA
  # skip this line if you run RStudio or other third-party R environments
  enableWGCNAThreads()
  
  # load data
  ds <- load.dataset(
    meta.file = './data/sample_sheet.csv', meta.sep = ',',
    data.file = './data/rna_norm_counts.csv', data.sep = ','
  )
  data <- ds$data.matrix
  
  n_samples <- nrow(data)
  
  # create topological overlap matrix
  if (file.exists('./rdata/wgcna.RData')) {
    print('loading wgcna data...')
    load('./rdata/wgcna.RData')
  } else {
    # search soft-thresholding powers
    powers = 2:20
    sft = pickSoftThreshold(data, powerVector = powers, verbose = 5)
    
    # plot the results
    sizeGrWindow(9, 5)
    par(mfrow = c(1,2))
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab = 'Soft Threshold (power)', ylab = 'Scale Free Topology Model Fit,signed R^2', type = 'n',
         main = paste('Scale independence'));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels = powers, cex = 0.9, col = 'red');
    abline(h = 0.90, col = 'red')
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab = 'Soft Threshold (power)', ylab = 'Mean Connectivity', type = 'n',
         main = paste('Mean connectivity'))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.9, col = 'red')
    
    # apply soft-thresholding
    softPower = 4 # chosen with elbow-method
    adjacency = adjacency(data, power = softPower)
    
    # topological similarity
    TOM = TOMsimilarity(adjacency)
    remove(adjacency)
    gc()
    
    dissTOM = 1 - TOM
    remove(TOM)
    gc()
    
    # hierarchical clustering for community detection
    geneTree = hclust(as.dist(dissTOM), method = 'average')
    
    save(dissTOM, geneTree, file = './rdata/wgcna.RData')
  }
  
  print('hierarchical clustering for community detection...')
  dynamicMods = cutreeDynamic(
    dendro = geneTree,
    distM = dissTOM,
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    minClusterSize = 20
  )
  moduleColors = labels2colors(dynamicMods)
  remove(dissTOM)
  gc()
  
  png(file=paste('results', 'wgcna', 'dendogram.png', sep = '/'), width = 1280)
  plotDendroAndColors(geneTree, moduleColors, 'Dynamic Tree Cut',
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  dev.off()
  
  print('calculate eigengenes...')
  MEList = moduleEigengenes(data, colors = moduleColors)
  MEs = MEList$eigengenes
  
  print('get conditions/responses...')
  conds <- lapply(ds$meta.data, function(x) as.numeric(as.factor(x)))
  conds <- as.data.frame(conds, row.names = row.names(ds$meta.data))
  
  print('done.')
  
  return(list(
    module_eigengenes = MEs,
    ds = ds,
    module_colors = moduleColors,
    conds = conds
  ))
}

cluster_module_eigengenes <- function (wgcna_result) {
  MEDiss = 1 - cor(wgcna_result$module_eigengenes);
  METree = hclust(as.dist(MEDiss), method = "average");
  
  png(file=paste('results', 'wgcna', 'eigengenes_clustering.png', sep = '/'), width = 860)
  plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
  dev.off()
}

significant_modules <- function (wgcna_result, response) {
  des_mat <- model.matrix(~ wgcna_result$ds$meta.data[[response]])
  fit <- limma::lmFit(t(wgcna_result$module_eigengenes), design = des_mat)
  fit <- limma::eBayes(fit)
  stats_df <- limma::topTable(fit, number = ncol(wgcna_result$module_eigengenes)) %>% tibble::rownames_to_column("module")
  
  most_significant_module <- stats_df[1, 'module']
  
  modules_df <- cbind(wgcna_result$module_eigengenes, wgcna_result$ds$meta.data)
  
  fname <- paste(paste('module_diff', tolower(response), sep = '_'), 'png', sep = '.')
  fname <- paste('results', 'wgcna', fname, sep = '/')
  plt <- ggplot(
    modules_df,
    aes(
      x = .data[[response]],
      y = modules_df[, most_significant_module],
      color = .data[[response]]
    )
  ) +
    ylab(most_significant_module) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    ggforce::geom_sina(maxwidth = 0.3) +
    theme_classic()
  ggsave(filename = fname, plot = plt, width = 8, height = 10)
  
  return(stats_df)
}

correlate_eigengenes_with_responses <- function(wgcna_result) {
  moduleTraitCor = cor(wgcna_result$module_eigengenes, wgcna_result$conds, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(wgcna_result$ds$data.matrix))
  
  # plot the results
  png(file=paste('results', 'wgcna', 'module_response_rel.png', sep = '/'), width = 1280, height = 720)
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(wgcna_result$conds),
                 yLabels = names(wgcna_result$module_eigengenes),
                 ySymbols = names(wgcna_result$module_eigengenes),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-response relationships"))
  dev.off()
}

get_module_membership <- function (wgcna_result) {
  data <- wgcna_result$ds$data.matrix
  n_samples <- nrow(data)
  module_names <- substring(names(wgcna_result$module_eigengenes), 3)
  
  gene_module_membership <- as.data.frame(cor(data, wgcna_result$module_eigengenes, use = "p"))
  names(gene_module_membership) <- paste("MM", module_names, sep="")
  
  mm_pvalue <- as.data.frame(corPvalueStudent(as.matrix(gene_module_membership), n_samples))
  names(mm_pvalue) <- paste("p.MM", module_names, sep="")
  
  gene_module_membership$gene <- rownames(gene_module_membership)
  rownames(gene_module_membership) <- NULL
  gene_module_membership <- gene_module_membership %>% dplyr::select(gene, everything())
  
  mm_pvalue$gene <- rownames(mm_pvalue)
  rownames(mm_pvalue) <- NULL
  mm_pvalue <- mm_pvalue %>% dplyr::select(gene, everything())
  
  return(list(
    gene_module_membership = gene_module_membership,
    mm_pvalue = mm_pvalue
  ))
}

get_gene_response_significance <- function (wgcna_result, response_str) {
  data <- wgcna_result$ds$data.matrix
  n_samples <- nrow(data)
  response = as.data.frame(wgcna_result$conds[[response_str]])
  names(response) <- response_str
  
  gene_response_significance <- as.data.frame(cor(data, response, use = "p"))
  names(gene_response_significance) <- paste("GS.", names(response), sep="")
  
  gs_pvalue <- as.data.frame(corPvalueStudent(as.matrix(gene_response_significance), n_samples))
  names(gs_pvalue) <- paste("p.GS.", names(response), sep="")
  
  gene_response_significance$gene <- rownames(gene_response_significance)
  rownames(gene_response_significance) <- NULL
  gene_response_significance <- gene_response_significance %>% dplyr::select(gene, everything())
  
  gs_pvalue$gene <- rownames(gs_pvalue)
  rownames(gs_pvalue) <- NULL
  gs_pvalue <- gs_pvalue %>% dplyr::select(gene, everything())
  
  return(list(
    gene_response_significance = gene_response_significance,
    gs_pvalue = gs_pvalue
  ))
}

get_mm_and_gs <- function (wgcna_result) {
  return(list(
    module_membership = get_module_membership(wgcna_result),
    gene_site_significance = get_gene_response_significance(wgcna_result, 'Site'),
    gene_ref_significance = get_gene_response_significance(wgcna_result, 'REF')
  ))
}

identify_interesting_genes <- function (wgcna_result, mm, gs, module) {
  module_names <- substring(names(wgcna_result$module_eigengenes), 3)
  column <- match(module, module_names)
  module_genes = wgcna_result$module_colors == module
  resp <- substring(colnames(gs$gene_response_significance)[2], 4)

  filename <- paste(paste('mm', 'gs', resp, module, sep = '_'), 'png', sep = '.')
  png(file=paste('results', 'wgcna', filename, sep = '/'), width = 1280, height = 720)
  par(mfrow = c(1,1))
  verboseScatterplot(abs(mm$gene_module_membership[module_genes, column + 1]),
                     abs(gs$gene_response_significance[module_genes, 2]),
                     xlab = paste('Module membership in', module, 'module'),
                     ylab = paste('Gene significance for', resp),
                     main = paste('Module membership vs. gene significance\n'),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
  dev.off()
  
  significant_genes_mm <- mm$mm_pvalue[module_genes,] %>%
    dplyr::select(gene, ends_with(paste('MM', module, sep = ''))) %>%
    arrange(pick(contains(module))) %>%
    filter(.[[2]] < 0.05)
  
  significant_genes_gs <- gs$gs_pvalue[module_genes,] %>%
    dplyr::select(gene, contains(resp)) %>%
    arrange(pick(contains(resp))) %>%
    filter(.[[2]] < 0.05)
  
  return(list(
    significant_genes_mm = significant_genes_mm,
    significant_genes_gs = significant_genes_gs
  ))
}

get_interesting_modules_results <- function (wgcna_result) {
  interesting_modules <- c('green', 'greenyellow', 'lightcyan', 'magenta', 'red', 'yellow', 'tan', 'cyan')
  interesting_modules_results <- list()
  
  mm_and_gs <- get_mm_and_gs(wgcna_result)
  
  for (mod in interesting_modules) {
    interesting_modules_results[[mod]] = list()
    
    interesting_genes_site <- identify_interesting_genes(
      wgcna_result,
      mm_and_gs$module_membership,
      mm_and_gs$gene_site_significance,
      mod
    )
    
    interesting_genes_site <- merge(
      interesting_genes_site$significant_genes_mm,
      interesting_genes_site$significant_genes_gs,
      by = 'gene'
    )
    
    interesting_genes_site_annotated <- annotate_genes(interesting_genes_site$gene, ont = 'BP')
    
    interesting_modules_results[[mod]][['interesting_genes_site']] = interesting_genes_site
    interesting_modules_results[[mod]][['interesting_genes_site_annotated']] = interesting_genes_site_annotated
    
    filename <- paste(paste('interesting_genes_site_annotated', mod, sep = '_'), 'png', sep = '.')
    path <- paste('results', 'wgcna', filename, sep = '/')
    plt <- goplot(interesting_genes_site_annotated)
    ggsave(filename = path, plot = plt, width = 10, height = 10)
    
    interesting_genes_ref <- identify_interesting_genes(
      wgcna_result,
      mm_and_gs$module_membership,
      mm_and_gs$gene_ref_significance,
      mod
    )
    
    interesting_genes_ref <- merge(
      interesting_genes_ref$significant_genes_mm,
      interesting_genes_ref$significant_genes_gs,
      by = 'gene'
    )
    
    interesting_genes_ref_annotated <- annotate_genes(interesting_genes_ref$gene, ont = 'BP')
    
    interesting_modules_results[[mod]][['interesting_genes_ref']] = interesting_genes_ref
    interesting_modules_results[[mod]][['interesting_genes_ref_annotated']] = interesting_genes_ref_annotated
    
    filename <- paste(paste('interesting_genes_ref_annotated', mod, sep = '_'), 'png', sep = '.')
    path <- paste('results', 'wgcna', filename, sep = '/')
    plt <- goplot(interesting_genes_ref_annotated)
    ggsave(filename = path, plot = plt, width = 10, height = 10)
  }
  
  return(interesting_modules_results)
}

# Export to cytoscape
# -------------------
# modules = c('cyan', 'tan') # select to modules of interest
# probes = colnames(data)
# inModule = is.finite(match(moduleColors, modules))
# modProbes = probes[inModule]
# modTOM = 1 - dissTOM[inModule, inModule]
# dimnames(modTOM) = list(modProbes, modProbes)
# efname = paste('./data/', paste(modules, collapse = '-'), '_edges.txt', sep = '')
# nfname = paste('./data/', paste(modules, collapse = '-'), '_nodes.txt', sep = '')
# 
# exportNetworkToCytoscape(modTOM,
#                          edgeFile = efname,
#                          nodeFile = nfname,
#                          weighted = TRUE, threshold = 0.25, 
#                          nodeNames = modProbes,
#                          nodeAttr = moduleColors[inModule])
