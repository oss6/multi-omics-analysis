library(clusterProfiler)
library(org.Dm.eg.db)

source('./id_mapping.R')

annotate_genes <- function (genes, ont) {
  id_map <- create_map(
    './annotations/rna_dma_to_dme_mappings.tsv',
    genes,
    all_mappings = TRUE
  )
  
  enrich_result <- enrichGO(id_map, 'org.Dm.eg.db', ont = ont, keyType = 'ENSEMBL', readable = TRUE)
  
  return(enrich_result)
}

annotate_genes_kegg <- function (genes) {
  ids <- create_map(
    './annotations/rna_dma_to_dme_mappings.tsv',
    genes,
    all_mappings = TRUE
  )
  
  ids_map <- bitr(ids, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Dm.eg.db')
  
  enrich_result <- enrichKEGG(ids_map$ENTREZID, organism = 'dme', keyType='ncbi-geneid')
  
  return(enrich_result)
}

clusters_comparison <- function (genes_list, kegg = TRUE) {
  gene_cluster <- lapply(genes_list, function (x) {
    ids <- create_map(
      './annotations/rna_dma_to_dme_mappings.tsv',
      x,
      all_mappings = TRUE
    )
    
    ids_map <- bitr(ids, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Dm.eg.db')
    
    return(ids_map$ENTREZID)
  })
  
  ck <- NULL
  
  if (kegg) {
    ck <- compareCluster(geneCluster = gene_cluster, fun = enrichKEGG, organism = 'dme', keyType = 'ncbi-geneid')
  } else {
    ck <- compareCluster(geneCluster = gene_cluster, fun = enrichGO, OrgDb = 'org.Dm.eg.db', keyType = 'ENTREZID')
  }
  
  ck <- setReadable(ck, OrgDb = 'org.Dm.eg.db', keyType = 'ENTREZID')
  
  return(ck)
}
