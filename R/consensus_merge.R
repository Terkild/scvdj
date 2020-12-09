#' Merge consensus VDJ annotation into contigs
#'
#' Merges consensus VDJ annotations into contig annotation data.frame to allow incomplete
#' to be called as most likely clonotype
#'
#' @param data_annotation Contig annotations from Cell Ranger output
#' @param data_consensus  Consensus contigs from Cell Ranger output
#' @param mutual_columns  Columns that are included from both outputs (set by mutual_columns)
#' @param consensus_suffix  Suffix given to mutual_columns from the consensus data
#'
#' @importFrom dplyr filter select inner_join one_of bind_rows
#' @return data.frame containing contig annotations
#' @export

consensus_merge <- function(data_annotation, data_consensus, mutual_columns=c("length","reads","umis"), consensus_suffix=".consensus"){
	## Determine which columns can be inherited from consensus annotation
	use.consensusColumns <- setdiff(intersect(colnames(data_consensus),colnames(data_annotation)),mutual_columns)

	data_annotation$raw_consensus_id <- gsub("consensus([0-9]+)$", "consensus_\\1", data_annotation$raw_consensus_id)

	## Replace contig annotation with consensus annotation for contigs assigned to a consensus
	vdj.consensus <- data_annotation %>%
		filter(raw_consensus_id != "None") %>%
		select(-one_of(use.consensusColumns)) %>%
		inner_join(data_consensus,
				   by=c("raw_consensus_id"="consensus_id"),
				   suffix=c("",consensus_suffix))

	## Get annotations lacking consensus
	vdj.noConsensus <- data_annotation %>% filter(raw_consensus_id == "None")

	## Recombine annotations
	vdj <- bind_rows(list(vdj.consensus,vdj.noConsensus))

	return(vdj)
}
