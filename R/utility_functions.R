## UTILITY FUNCTION

#' Determines if a field from cellranger vdj output is empty
#'
#' Considered empty if equal to "" or "None"
#'

no_empty <- function(x, empty.values=c("","None")){x[!x %in% empty.values]}


#' Make VDJ string
#'
#' Paste values of multiple columns into a single string
paste_vdj <- function(x, sep=":"){
	paste(apply(x,1,function(y){
		paste(no_empty(y),collapse=sep)
	}),collapse="; ")
}
