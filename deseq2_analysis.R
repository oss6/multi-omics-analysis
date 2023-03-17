library(DESeq2)
library(tidyverse)

source('./fileio.R')
source('./genes_annotation.R')

to_tibble <- function (.data) {
  data.frame(.data) %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
}

get_degs <- function () {
  rna_raw.ds <- load.dataset(
    meta.file = './data/sample_sheet.csv', meta.sep = ',',
    data.file = './data/rna_raw_counts.csv', data.sep = ','
  )
  
  rna_raw.ds$data.matrix <- round(rna_raw.ds$data.matrix)
  
  control_samples_removal <- -seq(1, 6)
  
  raw_counts <- t(rna_raw.ds$data.matrix)
  raw_counts <- raw_counts[,control_samples_removal]
  
  sample_sheet <- rna_raw.ds$meta.data
  sample_sheet <- sample_sheet[control_samples_removal,]
  sample_sheet$REF <- relevel(factor(sample_sheet$REF), ref = "1x")
  sample_sheet$Site <- relevel(factor(sample_sheet$Site), ref = "D01")
  
  # Remove low/zero variance genes
  read_vars <- rowVars(raw_counts)
  raw_counts <- raw_counts[which(read_vars > 1e-10),]
  
  # Remove genes with too many null counts
  read_valids <- rowSums(raw_counts != 0)
  raw_counts <- raw_counts[which(read_valids > 3),]
  
  # Remove genes with too small total counts
  read_cnts <- rowSums(raw_counts)
  raw_counts <- raw_counts[which(read_cnts > 200),]
  
  # Differential analysis
  dds <- DESeqDataSetFromMatrix(raw_counts, colData = sample_sheet, design = ~REF + Site + REF:Site)
  dds <- DESeq(dds)
  
  resultsNames(dds)
  
  res <- results(dds, name = 'REF_10x_vs_1x')
  
  degs <- to_tibble(res) %>% filter(padj < 0.05) %>% filter(abs(log2FoldChange) > 1)
  
  degs_annotated <- annotate_genes(degs$gene, ont = 'BP')
  
  return(list(
    degs = degs,
    degs_annotated = degs_annotated
  ))
}
