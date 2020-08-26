#' Add VDJ info to Seurat object
#'
#' Adds parsed contig by cell dataframe into metadata of a Seurat object
#'
#' @param object  Seurat object
#' @param VDJ Parsed VDJ contig by cell data.frame
#' @param column.CB Name of column that contains cell barcode
#' @param prefix  Prefix given to each metadata column
#' @param colums.include  Columns from parsed contig by cell data.frame to be included in Seurat metadata
#' @param chains  Vector of chains (TRA, TRB, TRG, TRD, IGH, IHK, IGL) to include from parsed contig by cell data.frame
#' @param create.boolean  Should a boolean of whether a given cell has any expression of a given chain be included
#' @param boolean.suffix  Suffix for boolean (default structure is CHAIN.bool)
#' @param boolean.filterFunctional Should boolean only consider "functional" VDJ counts
#'
#' @return Seurat object
#' @importFrom Seurat AddMetaData
#' @export

seurat_vdj_add_to_meta <- function(object, VDJ, column.CB="CB", prefix="VDJ_", columns.include=c("clonotype.1","clonotype.2","clonotype.other","CDR3.1","CDR3.2","CDR3.other","V.1","V.2","V.other","D.1","D.2","D.other","J.1","J.2","J.other","C.1","C.2","C.other","UMIcount.1","UMIcount.2","UMIcount.other","reads.1","reads.2","reads.other"), chains=NULL, create.boolean=TRUE, boolean.suffix=".bool", boolean.filterFunctional=FALSE){

	VDJ.byChain <- split(VDJ, f=VDJ$topChains)
	VDJ.byChain <- lapply(VDJ.byChain,function(x){rownames(x) <- x[,column.CB]; return(x)})
	VDJ.byChain <- lapply(VDJ.byChain,function(x)x[intersect(rownames(x),colnames(object)),])

	if(is.null(chains)){
		chains <- names(VDJ.byChain)
	}

	for(i in seq_along(chains)){
		curChain <- chains[i]

		for(j in seq_along(columns.include)){
			curColumn <- columns.include[j]

			col.name <- paste0(prefix,curChain,".",curColumn)

			object <- AddMetaData(object, NA, col.name=col.name)
			if(length(na.omit(VDJ.byChain[[curChain]][,curColumn]))==length(grep("^[0-9]+$",na.omit(VDJ.byChain[[curChain]][,curColumn])))){
				VDJ.byChain[[curChain]][,curColumn] <- as.numeric(VDJ.byChain[[curChain]][,curColumn])
				VDJ.byChain[[curChain]][is.na(VDJ.byChain[[curChain]][,curColumn]),curColumn] <- 0
			}
			object[[col.name]][VDJ.byChain[[curChain]][,column.CB],] <- VDJ.byChain[[curChain]][,curColumn]
		}

		if(create.boolean == TRUE){
			col.name <- paste0(prefix,curChain,boolean.suffix)
			object[[col.name]] <- FALSE

			if(boolean.filterFunctional == TRUE){
				## You can filter cells with nonfunctional CDR3
				object[[col.name]][VDJ.byChain[[curChain]][grep("\\*|_",VDJ.byChain[[curChain]]$CDR3.1,invert=TRUE),column.CB],] <- TRUE
			} else {
				object[[col.name]][VDJ.byChain[[curChain]][,column.CB],] <- TRUE
			}
		}
	}

	return(object)
}
