if (!("clusterProfiler" %in% installed.packages())) {
  BiocManager::install("clusterProfiler")
}

if (!("org.Dm.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Dm.eg.db")
}

if (!("FELLA" %in% installed.packages())) {
  BiocManager::install("FELLA")
}

install.packages("ggcorrplot")

install.packages("FactoMineR")

install.packages("factoextra")

remove.packages("RGCCA")
devtools::install_github(repo="https://github.com/rgcca-factory/RGCCA.git", ref = "main")

library(dplyr)
library(enrichplot)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)

if (!dir.exists('./rdata')) {
  dir.create('./rdata')
}

if (!dir.exists('./results')) {
  dir.create('./results')
}

if (!dir.exists('./results/rgcca')) {
  dir.create('./results/rgcca')
}

if (!dir.exists('./results/wgcna')) {
  dir.create('./results/wgcna')
}

if (!dir.exists('./results/irf')) {
  dir.create('./results/irf')
}

source('./compounds_annotation.R')
source('./genes_annotation.R')
source('./id_mapping.R')
source('./rgcca.R')
source('./irf.R')
source('./deseq2_analysis.R')
source('./wgcna.R')

fella_data <- load_compounds_annotation_data('dmk')

genes_characterising_site <- list()
genes_characterising_ref <- list()

compounds_characterising_site <- list()
compounds_characterising_ref <- list()

# DESeq2
# ------

degs_result <- get_degs()

genes_characterising_ref[['degs']] <- degs_result$degs$gene

# goplot(degs_result$degs_annotated)
# dotplot(degs_result$degs_annotated)

# R/SGCCA
# ------------------------------------------------------------------------------

# REF

ref_allpeaks_res <- run_rgcca(degs_result$degs$gene, 'ref', use_specific_peaks = FALSE)
ref_allpeaks_res_significant_genes <- ref_allpeaks_res$boot_res$stats %>%
  filter(th_pval < 0.05) %>%
  dplyr::select(var, block, comp, type) %>%
  filter(block == 'rna' & comp == 1 & type == 'weights') %>%
  distinct(var)
ref_allpeaks_res_significant_compounds <- ref_allpeaks_res$boot_res$stats %>%
  filter(th_pval < 0.05) %>%
  dplyr::select(var, block, comp, type) %>%
  filter(block == 'ppolar' & comp == 1 & type == 'weights') %>%
  distinct(var)

ref_specpeaks_res <- run_rgcca(degs_result$degs$gene, 'ref', use_specific_peaks = TRUE)
ref_specpeaks_res_significant_genes <- ref_specpeaks_res$boot_res$stats %>%
  filter(th_pval < 0.05) %>%
  dplyr::select(var, block, comp, type) %>%
  filter(block == 'rna' & comp == 1 & type == 'weights') %>%
  distinct(var)
ref_specpeaks_res_significant_compounds <- ref_allpeaks_res$boot_res$stats %>%
  filter(th_pval < 0.05) %>%
  dplyr::select(var, block, comp, type) %>%
  filter(block != 'rna' & comp == 1 & type == 'weights') %>%
  distinct(var)

genes_characterising_ref[['sgcca']] <- unique(c(ref_allpeaks_res_significant_genes$var, ref_specpeaks_res_significant_genes$var))
compounds_characterising_ref[['sgcca']] <- unique(c(ref_allpeaks_res_significant_compounds$var, ref_specpeaks_res_significant_compounds$var))

## Annotation

sgcca_genes_ref <- annotate_genes(genes_characterising_ref$sgcca, ont = 'BP')
p <- dotplot(sgcca_genes_ref, title = 'BP characterising REF')
ggsave(filename = './results/rgcca/genes_ref_bp_dotplot.png', plot = p, width = 6, height = 10)

sgcca_compounds_ref <- annotate_compounds(fella_data, compounds_characterising_ref$sgcca)
plot_compounds(sgcca_compounds_ref, './results/rgcca/significant_compounds_graph_ref.png')
save_compounds_table(sgcca_compounds_ref, './results/rgcca/significant_compounds_table_ref.csv')


# Site

site_allpeaks_res <- run_rgcca(degs_result$degs$gene, 'site', use_specific_peaks = FALSE)
site_allpeaks_res_significant_genes <- site_allpeaks_res$boot_res$stats %>%
  filter(th_pval < 0.05) %>%
  dplyr::select(var, block, comp, type) %>%
  filter(block == 'rna' & comp == 1 & type == 'weights') %>%
  distinct(var)
site_allpeaks_res_significant_compounds <- site_allpeaks_res$boot_res$stats %>%
  filter(th_pval < 0.05) %>%
  dplyr::select(var, block, comp, type) %>%
  filter(block != 'rna' & comp == 1 & type == 'weights') %>%
  distinct(var)

site_specpeaks_res <- run_rgcca(degs_result$degs$gene, 'site', use_specific_peaks = TRUE)
site_specpeaks_res_significant_genes <- site_specpeaks_res$boot_res$stats %>%
  filter(th_pval < 0.05) %>%
  dplyr::select(var, block, comp, type) %>%
  filter(block == 'rna' & comp == 1 & type == 'weights') %>%
  distinct(var)
site_specpeaks_res_significant_compounds <- site_specpeaks_res$boot_res$stats %>%
  filter(th_pval < 0.05) %>%
  dplyr::select(var, block, comp, type) %>%
  filter(block != 'rna' & comp == 1 & type == 'weights') %>%
  distinct(var)

genes_characterising_site[['sgcca']] <- unique(c(site_allpeaks_res_significant_genes$var, site_specpeaks_res_significant_genes$var))
compounds_characterising_site[['sgcca']] <- unique(c(site_allpeaks_res_significant_compounds$var, site_specpeaks_res_significant_compounds$var))

## Annotation

sgcca_genes_site <- annotate_genes(genes_characterising_site$sgcca, ont = 'BP')
p <- dotplot(sgcca_genes_site, title = 'BP characterising Site')
ggsave(filename = './results/rgcca/genes_site_bp_dotplot.png', plot = p, width = 6, height = 10)

sgcca_compounds_site <- annotate_compounds(fella_data, compounds_characterising_site$sgcca)
plot_compounds(sgcca_compounds_site, './results/rgcca/significant_compounds_graph_site.png')
save_compounds_table(sgcca_compounds_site, './results/rgcca/significant_compounds_table_site.csv')

# ------------------------------------------------------------------------------



# WGCNA
# ------------------------------------------------------------------------------

wgcna_result <- run_wgcna()
cluster_module_eigengenes(wgcna_result)
site_significant_modules <- significant_modules(wgcna_result, 'Site')
ref_significant_modules <- significant_modules(wgcna_result, 'REF')
correlate_eigengenes_with_responses(wgcna_result)
interesting_modules <- get_interesting_modules_results(wgcna_result)

wgcna_site_genes <- c()
wgcna_ref_genes <- c()

for (mod in interesting_modules) {
  wgcna_site_genes <- c(wgcna_site_genes, mod$interesting_genes_site$gene)
  wgcna_ref_genes <- c(wgcna_ref_genes, mod$interesting_genes_ref$gene)
}

genes_characterising_site[['wgcna']] <- wgcna_site_genes
genes_characterising_ref[['wgcna']] <- wgcna_ref_genes

wgcna_genes_site_bp <- annotate_genes(genes_characterising_site$wgcna, ont = 'BP')
ggsave(
  filename = './results/wgcna/genes_site_bp_dotplot.png',
  plot = dotplot(wgcna_genes_site_bp, title = 'BP characterising Site'),
  width = 6,
  height = 10)

wgcna_genes_site_mf <- annotate_genes(genes_characterising_site$wgcna, ont = 'MF')
ggsave(
  filename = './results/wgcna/genes_site_mf_dotplot.png',
  plot = dotplot(wgcna_genes_site_mf, title = 'MF characterising Site'),
  width = 6,
  height = 10)

wgcna_genes_ref_bp <- annotate_genes(genes_characterising_ref$wgcna, ont = 'BP')
ggsave(
  filename = './results/wgcna/genes_ref_bp_dotplot.png',
  plot = dotplot(wgcna_genes_ref_bp, title = 'BP characterising REF'),
  width = 6,
  height = 10)

wgcna_genes_ref_mf <- annotate_genes(genes_characterising_ref$wgcna, ont = 'MF')
ggsave(
  filename = './results/wgcna/genes_ref_mf_dotplot.png',
  plot = dotplot(wgcna_genes_ref_mf, title = 'MF characterising REF'),
  width = 6,
  height = 10)

# ------------------------------------------------------------------------------



# iRF
# ------------------------------------------------------------------------------

# Get significant compounds against REF

compounds_ref_irf_result <- get_important_features.irf('./data/polar_pos_pqn_imputed_glog.csv', 'REF', ref_level = '1x')
plot_stability_scores(compounds_ref_irf_result$irf_fit, './results/irf/ref_compound_interactions_stability_scores.png')
compounds_ref_irf_significant_peaks <- rownames(compounds_ref_irf_result$significant_features)
compounds_ref_irf_analysis <- annotate_compounds(fella_data, compounds_ref_irf_significant_peaks)
plot_compounds(compounds_ref_irf_analysis, './results/irf/ref_significant_compounds.png')
save_compounds_table(compounds_ref_irf_analysis, './results/irf/ref_significant_compounds.csv')

compounds_characterising_ref[['irf']] <- compounds_ref_irf_significant_peaks


# Get significant compounds against Site

compounds_site_irf_result <- get_important_features.irf('./data/polar_pos_pqn_imputed_glog.csv', 'Site', ref_level = 'D01')
plot_stability_scores(compounds_site_irf_result$irf_fit, './results/irf/site_compound_interactions_stability_scores.png')
compounds_site_irf_significant_peaks <- rownames(compounds_site_irf_result$significant_features)
compounds_site_irf_analysis <- annotate_compounds(fella_data, compounds_site_irf_significant_peaks)
plot_compounds(compounds_site_irf_analysis, './results/irf/site_significant_compounds.png')
save_compounds_table(compounds_site_irf_analysis, './results/irf/site_significant_compounds.csv')

compounds_characterising_site[['irf']] <- compounds_site_irf_significant_peaks

# Get significant genes against REF

genes_ref_irf_result <- get_important_features.irf('./data/rna_norm_counts.csv', 'REF', ref_level = '1x', mask = degs_result$degs$gene)
plot_stability_scores(genes_ref_irf_result$irf_fit, './results/irf/ref_genes_interactions_stability_scores.png')

genes_characterising_ref[['irf']] <- rownames(genes_ref_irf_result$significant_features)

irf_genes_ref <- annotate_genes(genes_characterising_ref$irf, ont = 'BP')
p <- dotplot(irf_genes_ref, title = 'BP characterising REF')
ggsave(filename = './results/irf/genes_ref_bp_dotplot.png', plot = p, width = 6, height = 10)

# Get significant genes against Site

genes_site_irf_result <- get_important_features.irf('./data/rna_norm_counts.csv', 'Site', ref_level = 'D01', mask = degs_result$degs$gene)
plot_stability_scores(genes_site_irf_result$irf_fit, './results/irf/site_genes_interactions_stability_scores.png')

genes_characterising_site[['irf']] <- rownames(genes_site_irf_result$significant_features)

irf_genes_site <- annotate_genes(genes_characterising_site$irf, ont = 'BP')
p <- dotplot(irf_genes_site, title = 'BP characterising Site')
ggsave(filename = './results/irf/genes_site_bp_dotplot.png', plot = p, width = 6, height = 10)



# ------------------------------------------------------------------------------


# Consensus
# ------------------------------------------------------------------------------

# Genes characterising Site

site_consensus_genes <- Reduce(union, genes_characterising_site)
site_consensus_genes_annotation_bp <- annotate_genes(site_consensus_genes, 'BP')
site_consensus_genes_annotation_mf <- annotate_genes(site_consensus_genes, 'MF')

ggsave(
  filename = './results/genes_site_bp_dotplot.png',
  plot = dotplot(site_consensus_genes_annotation_bp, title = 'BP characterising Site'),
  width = 6,
  height = 10)

ggsave(
  filename = './results/genes_site_mf_dotplot.png',
  plot = dotplot(site_consensus_genes_annotation_mf, title = 'MF characterising Site'),
  width = 6,
  height = 10)

# Genes characterising REF

ref_consensus_genes <- Reduce(union, genes_characterising_ref)
ref_consensus_genes_annotation_bp <- annotate_genes(ref_consensus_genes, 'BP')
ref_consensus_genes_annotation_mf <- annotate_genes(ref_consensus_genes, 'MF')

ggsave(
  filename = './results/genes_ref_bp_dotplot.png',
  plot = dotplot(ref_consensus_genes_annotation_bp, title = 'BP characterising REF'),
  width = 6,
  height = 10)

ggsave(
  filename = './results/genes_site_mf_dotplot.png',
  plot = dotplot(ref_consensus_genes_annotation_mf, title = 'MF characterising REF'),
  width = 6,
  height = 10)

# Gene comparison between Site and REF

source('./genes_annotation.R')

genes_list <- list(site = genes_characterising_site$sgcca, ref = genes_characterising_ref$sgcca)
ck <- clusters_comparison(genes_list, kegg = FALSE)
ck <- pairwise_termsim(ck)
svg(filename = './results/sgcca_site_vs_ref_emapplot.svg')
emapplot(ck)
dev.off()

ck <- clusters_comparison(genes_list, kegg = TRUE)
ck <- pairwise_termsim(ck)
svg(filename = './results/sgcca_site_vs_ref_kegg_emapplot.svg')
emapplot(ck)
dev.off()

# Gene comparison between methods

ck <- clusters_comparison(genes_characterising_site, kegg = FALSE)
ck <- pairwise_termsim(ck)
svg(filename = './results/site_methods_emapplot.svg', width = 15, height = 15)
emapplot(ck)
dev.off()

# Compounds characterising Site

site_consensus_compounds <- Reduce(union, compounds_characterising_site)
site_consensus_compounds_annotated <- annotate_compounds(fella_data, site_consensus_compounds, nlimit_t = 500)

plot_compounds(site_consensus_compounds_annotated, './results/compounds_site.png')
save_compounds_table(site_consensus_compounds_annotated, './results/compounds_site.csv')

# Compounds characterising REF

ref_consensus_compounds <- Reduce(union, compounds_characterising_ref)
ref_consensus_compounds_annotated <- annotate_compounds(fella_data, ref_consensus_compounds, nlimit_t = 500)

plot_compounds(ref_consensus_compounds_annotated, './results/compounds_ref.png')
save_compounds_table(ref_consensus_compounds_annotated, './results/compounds_ref.csv')


# Water chemicals analysis
# ------------------------------------------------------------------------------

water_chemicals <- read.table('./data/water_chemicals.tsv', header = T, sep = '\t', row.names = 1, check.names = F)
water_chemicals <- as.data.frame(water_chemicals[,-c(1)])
# water_chemicals_meta <- data.frame(name = rownames(water_chemicals), cas = water_chemicals$CAS)
water_chemicals <- as.data.frame(t(water_chemicals[, -1]))
water_chemicals <- water_chemicals %>% select_if(colSums(.) != 0)
water_chemicals <- scale(water_chemicals)

water_chemicals_corr <- cor(water_chemicals, use = 'pairwise.complete.obs')
water_chemicals_pca <- princomp(water_chemicals_corr)

fviz_eig(water_chemicals_pca, addlabels = TRUE)

# fviz_pca_var(water_chemicals_pca, col.var = "black")
# 
# fviz_cos2(water_chemicals_pca, choice = "var", axes = 1:2)

fviz_pca_var(water_chemicals_pca, col.var = "cos2",
             gradient.cols = c("red", "gray", "blue"),
             repel = TRUE)



# ---------------------


sites <- unlist(lapply(rep(1:12, each = 6), function (x) { paste('D', formatC(x, width = 2, flag = "0"), sep = '') }))

g <- as.data.frame(rna.ds$data.matrix) %>% select(any_of(site_consensus_genes))
g <- g[-seq(1, 6),]
g$Sample <- rownames(g)
g$Site <- c(sites, sites)

chem_df <- read.table('./data/water_chemicals.tsv', header = T, sep = '\t', row.names = 1, check.names = F)
chem_df <- chem_df[,-1]

# Convert the chemicals matrix to long format
melt_chem_df <- melt(chem_df, value.name = "Concentration") %>% rename(Chemical = ChemName, Site = variable)

# Merge the normalized gene expression data with the long format chemicals matrix
merged_df <- merge(g, melt_chem_df, by = 'Site')

# Create a design matrix for the linear model
#design_matrix <- model.matrix(~ 0 + site + factor(Sample_ID))

design_matrix <- model.matrix(~ Chemical + Site, merged_df)

merged_df <- merged_df %>% select(-one_of(c('Chemical', 'Site', 'Sample')))
merged_df$Site <- factor(merged_df$Site)
merged_df$Chemical <- factor(merged_df$Chemical)

# Fit the linear model
fit <- lmFit(merged_df %>% select(-one_of(c('Chemical', 'Site', 'Sample'))), design_matrix)

# Apply empirical Bayes smoothing to the model
fit <- eBayes(fit)

# Identify differentially expressed genes
results <- topTable(fit, coef = "concentration", number = Inf, adjust.method = "BH", sort.by = "B", p.value = 0.05)



# # Read in the three datasets
# water_chem <- read.csv("water_chem_data.csv")
# gene_expression <- read.csv("gene_expression_data.csv")
# metabolite_data <- read.csv("metabolite_data.csv")
# 
# # Reshape the water chemicals matrix into a long format
# library(reshape2)
# water_chem_long <- reshape2::melt(water_chemicals, id.vars = "Site", variable.name = "Chemical", value.name = "Concentration") %>% rename(Site = Var1, Chemical = Var2)
# 
# # Normalize the concentration values
# # reference_ion <- "Na+"  # choose a commonly found ion as the reference
# # water_chem_long$Concentration_norm <- water_chem_long$Concentration / water_chem_long[water_chem_long$Chemical == reference_ion, "Concentration"]
# 
# # Reshape the gene expression and metabolite matrices into a long format
# sites <- unlist(lapply(rep(1:12, each = 6), function (x) { paste('D', formatC(x, width = 2, flag = "0"), sep = '') }))
# rna.ds <- load.dataset(
#   meta.file = './data/sample_sheet.csv', meta.sep = ',',
#   data.file = './data/rna_norm_counts.csv', data.sep = ','
# )
# g <- as.data.frame(rna.ds$data.matrix) %>% select(any_of(genes_characterising_site$wgcna))
# g <- g[-seq(1, 6),]
# g$Sample <- rownames(g)
# g$Site <- c(sites, sites)
# gene_expression_long <- melt(g, id.vars = c("Sample", "Site"), variable.name = "Gene", value.name = "Expression")
# metabolite_data_long <- melt(metabolite_data, id.vars = "Sample", variable.name = "Metabolite", value.name = "Intensity")
# 
# # Merge the three datasets by sample ID
# library(dplyr)
# merged_data <- inner_join(gene_expression_long, metabolite_data_long, by = "Sample")
# merged_data <- inner_join(gene_expression_long, water_chem_long, by = "Site")
# 
# # Perform correlation analysis
# gene_corr <- cor(merged_data$Expression, merged_data$Concentration, method = "pearson")
# metabolite_corr <- cor(merged_data$Intensity, merged_data$Concentration_norm, method = "pearson")
# 
# # Visualize the results
# plot(merged_data$Concentration_norm, merged_data$Expression, xlab = "Normalized Concentration", ylab = "Gene Expression")
# plot(merged_data$Concentration_norm, merged_data$Intensity, xlab = "Normalized Concentration", ylab = "Metabolite Intensity")
