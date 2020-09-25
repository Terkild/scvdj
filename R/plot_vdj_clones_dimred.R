#' Plot VDJ clones on DimPlot seurat
#'
#' @param object Seurat object
#'
#' @return ggplot object
#' @export
#' @importFrom Seurat FetchData

plot_vdj_clones_dimred_seurat <- function(object, clone, reduction="umap", ...){
  red.key <- object@reductions[[reduction]]@key
  red.dims <- paste0(red.key,c(1,2))

  data <- Seurat::FetchData(object, c(clone, red.dims))

  colnames(data) <- c("clone", "dim1", "dim2")

  return(plot_vdj_clones_dimred(x=data$dim1, y=data$dim2, clone=data$clone, red.dims=red.dims...))
}


#' Plot VDJ clones on DimPlot
#'
#'
#' @import ggplot2

#' @return ggplot object
#' @export
#'

# Redo function to be Seurat independent and make wrapper function for seurat
plot_vdj_clones_dimred <- function(x, y, clone, clone_colors=c(), red.dims=c("dim1","dim2"), polyclonal_label="Polyclonal", polyclonal_color="#6666FF", clone_size=2, polyclonal_size=0.75){
  plotData <- data.frame(x=x, y=y, clone=clone)

  p <- ggplot(plotData, aes(x=x, y=y)) +
    geom_point(color="lightgrey", size=polyclonal_size) +
    geom_point(data=plotData[which(plotData$clone == polyclonal_label),], color=polyclonal_color, size=polyclonal_size) +
    geom_point(data=plotData[-which(plotData$clone %in% c(NA,polyclonal_label)),],
               pch=21, color="black", aes(fill=clone)) +
    guides(fill=guide_legend(ncol=1, override.aes=list(size=clone_size))) +
    labs(x=red.dims[1], y=red.dims[2]) +
    theme(legend.text=element_text(size=6),
          legend.key.height=unit(2,"mm"),
          plot.background=element_blank())

  if(length(clone_colors) >= length(unique(clone))) p <- p + scale_fill_manual(values=clone_colors)

  return(p)
}
