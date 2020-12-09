#' Load Cell Ranger VDJ
#'
#'
#' @param path	Path to Cell Ranger output (root folder, not "outs" folder)
#' @param consensus Boolean if contigs assigned to consensus sequences should use the consensus annotation
#' @param annotation_file contig_annotation file to use ("all" versus "filtered")
#' @param consensus_file  consensus_annotation file to use
#' @param ...	Passes variables on to \code{\link{consensus_merge}}
#'
#' @return data.frame containing contig annotations
#' @export

cellranger_vdj_load <- function(path, consensus=TRUE, annotation_file="outs/filtered_contig_annotations.csv", consensus_file="outs/consensus_annotations.csv", ...){

  data.annotation <- read.csv(file.path(path,annotation_file), sep=",", header=TRUE)

  if(consensus == TRUE){
    data.consensus <- read.csv(file.path(path,consensus_file), sep=",", header=TRUE)

    vdj <- consensus_merge(data.annotation, data.consensus, ...)

    ## Test if something wierd has happned
    #nrow(vdj) == nrow(data.annotation)
  } else {
    vdj <- data.annotation
    vdj$clonotype_id <- vdj$raw_clonotype_id
  }

  return(vdj)
}

#' Filter Cell Ranger VDJ
#'
#' Filter by boolean columns
#'
#' @param contig_annotations  Dataframe containing contig_annotations
#' @param filterTrue  Vector of which "Boolean columns" in cellranger contig_annotation should be true ("full_length", "high_confidence", "productive")
#' @param min.umis  Integer threshold for minimum UMI count required to be kept
#'
#' @importFrom dplyr filter_at all_of all_vars filter
#' @return data.frame containing contig annotations
#' @export

cellranger_vdj_filter <- function(contig_annotations, filterTrue=c("full_length","high_confidence"), min.umis=0){
  contig_annotations %>% filter_at(all_of(filterTrue), all_vars(. == "true")) %>% filter(umis >= min.umis)
}

#' Order Cell Ranger VDJ calls
#'
#' To determine which contigs are most likely to be "primary" (~expressed) loci
#'
#' @param contig_annotations  Dataframe containing contig_annotations
#'
#' @importFrom dplyr group_by arrange desc mutate
#' @return data.frame containing contig annotations
#' @export

cellranger_vdj_order <- function(contig_annotations){
  ## Make an order variable determined by being having functional rearranged CDR3, most UMIs or most reads
  alignments.byCell.byV.clonotypes <- contig_annotations %>%
    group_by(barcode, chain) %>%
    arrange(barcode, chain, desc(as.integer(umis)), clonotype_id, desc(as.integer(reads))) %>%
    mutate(ord=seq_along(cdr3), topChains=chain, CB=barcode)

  return(alignments.byCell.byV.clonotypes)
}

#' Parse cellranger contigs
#'
#' Parse cellranger contigs contig_annotation table to call VDJ on a per
#' cell-barcode basis
#'
#' @param contig_annotations  Dataframe containing contig_annotations
#'
#' @importFrom dplyr group_by summarize
#' @return data.frame containing contig annotations
#' @export

cellranger_vdj_parse_by_cell <- function(contig_annotations){

  ## There are some cells with multiple alleles of the same chain
  clonotype.alleles.byCell <- contig_annotations %>%
    mutate(UMIcount=as.character(umis), reads=as.character(reads)) %>%
    group_by(CB, topChains) %>%
    summarize(clonotype.1=clonotype_id[ord==1],
              clonotype.2=ifelse(max(ord)>1,clonotype_id[ord==2],""),
              clonotype.other=ifelse(max(ord)>2,paste(clonotype_id[ord>2], collapse="; "),""),
              CDR3.1=cdr3[ord==1],
              CDR3.2=ifelse(max(ord)>1,cdr3[ord==2],""),
              CDR3.other=ifelse(max(ord)>2,paste(cdr3[ord>2], collapse="; "),""),
              UMIcount.1=UMIcount[ord==1],
              UMIcount.2=ifelse(max(ord)>1,UMIcount[ord==2],""),
              UMIcount.other=ifelse(max(ord)>2,paste(UMIcount[ord>2], collapse="; "),""),
              reads.1=reads[ord==1],
              reads.2=ifelse(max(ord)>1,reads[ord==2],""),
              reads.other=ifelse(max(ord)>2,paste(reads[ord>2], collapse="; "),""),
              V.1=v_gene[ord==1],
              V.2=ifelse(max(ord)>1,v_gene[ord==2],""),
              V.other=ifelse(max(ord)>2,paste(v_gene[ord>2], collapse="; "),""),
              D.1=d_gene[ord==1],
              D.2=ifelse(max(ord)>1,d_gene[ord==2],""),
              D.other=ifelse(max(ord)>2,paste(d_gene[ord>2], collapse="; "),""),
              J.1=j_gene[ord==1],
              J.2=ifelse(max(ord)>1,j_gene[ord==2],""),
              J.other=ifelse(max(ord)>2,paste(j_gene[ord>2], collapse="; "),""),
              C.1=c_gene[ord==1],
              C.2=ifelse(max(ord)>1,c_gene[ord==2],""),
              C.other=ifelse(max(ord)>2,paste(c_gene[ord>2], collapse="; "),""))

  return(clonotype.alleles.byCell)
}
