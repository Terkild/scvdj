## UTILITY FUNCTION

#' Determines if a field from cellranger vdj output is empty
#'
#' Considered empty if equal to "" or "None"
#'
#' @param x Input vector
#' @param empty.values  Vector of values considered "empty"
#'
#' @return Vector where empty values have been removed
#' @export

no_empty <- function(x, empty.values=c("","None")){x[!x %in% empty.values]}


#' Make VDJ string
#'
#' Concatenate values of multiple columns into a single string after removing empty values by \code{no_empty}
#'
#' @param x vector of values
#' @param sep seperator in concatenated string
#'
#' @return Concatenated string where empty values have been removed
#' @export

paste_vdj <- function(x, sep=":"){
	paste(apply(x,1,function(y){
		paste(no_empty(y),collapse=sep)
	}),collapse="; ")
}
