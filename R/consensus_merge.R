#' Merge consensus VDJ annotation into contigs
#'
#' Merges consensus VDJ annotations into contig annotation data.frame to allow incomplete
#' to be called as most likely clonotype
#'

consensus_merge <- function(data.annotation, data.consensus, consensus.suffix=".consensus"){
	## Determine which columns can be inherited from consensus annotation
	use.consensusColumns <- setdiff(intersect(colnames(data.consensus),colnames(data.annotation)),c("length","reads","umis"))

	## Replace contig annotation with consensus annotation for contigs assigned to a consensus
	vdj.consensus <- data.annotation %>%
		filter(raw_consensus_id != "None") %>%
		select(-one_of(use.consensusColumns)) %>%
		inner_join(data.consensus,
				   by=c("raw_consensus_id"="consensus_id"),
				   suffix=c("",consensus.suffix))

	## Get annotations lacking consensus
	vdj.noConsensus <- data.annotation %>% filter(raw_consensus_id == "None")

	## Recombine annotations
	vdj <- bind_rows(list(vdj.consensus,vdj.noConsensus))

	return(vdj)
}
