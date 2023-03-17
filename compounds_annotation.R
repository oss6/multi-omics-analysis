library(FELLA)

load_compounds_annotation_data <- function (organism) {
  dbdir <- './keggdb'
  
  if (!dir.exists(dbdir)) {
    graph <- buildGraphFromKEGGREST(organism = organism)
    
    buildDataFromGraph(
      keggdata.graph = graph,
      databaseDir = dbdir,
      internalDir = FALSE,
      matrices = "diffusion",
      normality = "diffusion",
      niter = 50)
  }
  
  fella.data <- loadKEGGdata(
    databaseDir = dbdir,
    internalDir = FALSE,
    loadMatrix = "diffusion")
  
  return(fella.data)
}

annotate_compounds <- function (fella_data, compounds, nlimit_g = 150, nlimit_t = 200) {
  comps <- create_map(
    './annotations/polar_pos_pkl_to_kegg_annotations.tsv',
    compounds,
    all_mappings = TRUE
  )
  
  analysis <- defineCompounds(compounds = comps, data = fella_data)
  analysis <- runDiffusion(
    object = analysis,
    data = fella_data,
    approx = "normality")
  
  results_graph <- generateResultsGraph(
    object = analysis,
    method = "diffusion",
    nlimit = nlimit_g,
    data = fella_data)
  
  results_table <- generateResultsTable(
    method = "diffusion",
    nlimit = nlimit_t,
    object = analysis,
    data = fella_data)
  
  return(list(
    analysis = analysis,
    results_graph = results_graph,
    results_table = results_table
  ))
}

plot_compounds <- function (analysis, filename) {
  png(filename = filename, width = 1920, height = 1200)
  plotGraph(analysis$results_graph, vertex.label.cex = .9, width = 1920, height = 1200)
  dev.off()
}

save_compounds_table <- function (analysis, filename) {
  write.csv(analysis$results_table %>% arrange(p.score), file = filename)
}
