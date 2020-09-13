#cdr3_functions
#' Filter CDR3 clones
#'
#' To find the most abundant CDR3 sequences
#'
#'

cdr3_clone_filter <- function(cdr3, cells_min=5, nonclonal_name="Polyclonal"){
  clones <- cdr3 %>% group_by(cdr3) %>% summarize()
}
