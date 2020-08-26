#' Plot VDJ clones on DimPlot
#'
#' @return ggplot object
plot_vdj_clones_dimred <- function(object, chain="TCRab_TRB", reduction="tsne", colors=NULL){
  red.key <- object@reductions[[reduction]]@key
  red.dims <- paste0(red.key,c(1,2))

  object[[paste0(chain,'.clone')]] <- ifelse(as.vector(object[[paste0(chain,'.CDR3.1')]][,1]) %in% names(colors), as.vector(object[[paste0(chain,'.CDR3.1')]][,1]), "Polyclonal")

  object[[paste0(chain,'.clone')]][object[[paste0(chain,'.bool')]] == 0] <- NA

  plotData <- FetchData(object, vars=c(red.dims,paste0(chain,".clone"),"sampleName"))
  colnames(plotData)[1:3] <- c("dim1","dim2","clone")

  p <- ggplot(plotData, aes(x=dim1, y=dim2)) +
    geom_point(color="lightgrey", size=0.75) +
    geom_point(data=plotData[which(plotData$clone == "Polyclonal"),], color="#6666FF", size=0.75) +
    geom_point(data=plotData[-which(plotData$clone %in% c(NA,"Polyclonal")),],
               pch=21, color="black", aes(fill=clone)) +
    guides(fill=guide_legend(ncol=1, override.aes=list(size=2))) +
    labs(title=chain, x=red.dims[1], y=red.dims[2]) +
    theme(legend.text=element_text(size=6),
          legend.key.height=unit(2,"mm"),
          plot.background=element_blank())

  if(!is.null(colors)) p <- p + scale_fill_manual(values=colors)
  return(p)
}
