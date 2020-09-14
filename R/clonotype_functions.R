#clonotype_functions

#' Filter clonotypes by number of cells
#'
#' @param clonotypes  Vector of clonotype assignment to each cell
#' @param cells_min Integer of number of cells required to be considered a clone
#' @param nonclonal_name  Character string for labelling clonotypes that do not fulfill filtering criteria
#'
#' @return Returns a vector of same length as input clonotypes with filtered clonotype assignments
#' @import dplyr
#' @export

clonotype_filter <- function(clonotypes, cells_min=5, nonclonal_name="Polyclonal"){
  clones <- data.frame(clone=clonotypes) %>%
    group_by(clone) %>% mutate(num_cells=n()) %>% ungroup() %>%
    mutate(clone_assignment=ifelse(num_cells >= cells_min, clone, nonclonal_name))

  return(clones$clone_assignment)
}

#' Filter top clonotypes by rank
#'
#' @param clonotypes  Vector of clonotype assignment to each cell
#' @param clonotypes_max Max number of clonotypes to return
#' @param cells_min Integer of number of cells required to be considered a clone
#' @param ... Passed on to clonotype_filter

clonotype_top <- function(clonotypes, clonotypes_max=10, cells_min=1, nonclonal_name="Polyclonal", ...){
  clones <- data.frame(clonotypes=clonotype_filter(clonotypes, cells_min=cells_min, ...)) %>%
    filter(clonotypes != nonclonal_name) %>%
    group_by(clonotypes) %>% summarize(num_cells=n()) %>%
    top_n(clonotypes_max)

  clonotypes <- data.frame(clonotypes=clonotypes) %>%
    left_join(clones) %>%
    mutate(clone=ifelse(is.na(num_cells),nonclonal_name,clonotypes))

  return(clonotypes$clone)
}
