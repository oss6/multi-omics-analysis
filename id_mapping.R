library(dplyr)

create_map <- function (mappings_file, ids, all_mappings = FALSE, save_target = NA) {
  id.maps <- read.csv(mappings_file, header = F, row.names = 1, sep = '\t', stringsAsFactors = F)
  
  id.maps.matched <- subset(id.maps, row.names(id.maps) %in% ids)
  id.maps.matched$orig <- row.names(id.maps.matched)
  row.names(id.maps.matched) <- NULL
  
  if (ncol(id.maps.matched) == 3) {
    id.maps.matched <- id.maps.matched %>% dplyr::rename(target = V3)
  } else {
    id.maps.matched <- id.maps.matched %>% dplyr::rename(target = V2)
  }
  
  id.maps.matched <- id.maps.matched[id.maps.matched$target != '',]
  
  if (all_mappings) {
    result <- unique(unlist(strsplit(id.maps.matched$target, split = ';')))
    
    if (!is.na(save_target)) {
      writeLines(result, save_target)
    }
  } else {
    result <- id.maps.matched %>% pmap_dfr(function(...) {
      current <- tibble(...)
      current %>% mutate(target = unique(unlist(strsplit(current$target, split = ';')))[1])
    })
    
    if (!is.na(save_target)) {
      writeLines(result$target, save_target)
    }
  }
  
  return(result)
}
