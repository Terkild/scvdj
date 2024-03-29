#' Circos plot of V and J segments
#'
#'
#' @param Vsegment  Vector of V segment assignment
#' @param Jsegment  Vector of J segment assignment
#' @param Dsegment  Vector of D segment assignment
#' @param clone Vector of clone assignment
#' @param polyclonal_label  What label is used to tag polyclonal (within clone vector)
#' @param polyclonal_color  Color for polyclonal connections
#' @param start_degree  Where on the circle should the first segment start
#' @param clone_opacity Opacity of clonal connections
#' @param clone_colors  Named vector of colors (names corresponds to values in clone vector)
#' @param link_lwd  Line width of connections
#'
#' @export
#' @import circlize
#' @import dplyr

plot_vdj_circos <- function(Vsegment, Jsegment, Dsegment=c(), clone=c(), polyclonal_label="Polyclonal", polyclonal_color="#00000020", start_degree=270, clone_opacity=0.8, clone_colors, link_lwd=0.25, label_threshold=2, remove_chain_pattern="^(IG[LKH])|^(TR[ABGD])"){

  data <- data.frame(clone=clone,
                     V=gsub(remove_chain_pattern,"",Vsegment),
                     J=gsub(remove_chain_pattern,"",Jsegment))

  data <- data %>%
    group_by(clone, V, J) %>%
    summarize(Freq=n())

  data <- data %>% arrange(desc(Freq))
  data <- data %>% filter(!is.na(V))

  clone_colors2 <- alpha(clone_colors,clone_opacity)
  names(clone_colors2) <- names(clone_colors)
  clone_colors2[polyclonal_label] <- polyclonal_color
  col <- clone_colors2[data$clone]

  border <- ifelse(data$clone != polyclonal_label, 1, NA)

  grid.col <- append(unique(data$V) %>% sort() %>% setNames(viridis::cividis(length(.)), .),
                     unique(data$J) %>% sort() %>% setNames(viridis::magma(length(.)), .))

  Jsegment.order <- data %>% group_by(J) %>% summarize(count=sum(Freq)) %>% arrange(count)
  Vsegment.order <- data %>% group_by(V) %>% summarize(count=sum(Freq)) %>% arrange(count)

  circos.clear()
  circos.par(start.degree=start_degree)

  chordDiagram(data[,c(2:4)], col=col, grid.col=grid.col, link.border=border, link.lwd=link_lwd,
               annotationTrack="grid", preAllocateTracks=list(track.height=max(strwidth(unlist(dimnames(data))))),
               #order=union(Vsegment.order[[1]],Jsegment.order[[1]]),
               link.largest.ontop=TRUE)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xplot = get.cell.meta.data("xplot")

    if(abs(xplot[2] - xplot[1]) > label_threshold) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
    }
  }, bg.border = NA) # here set bg.border to NA is important

}

#' Circos plot of V and J segments from Seurat object
#'
#'
#' @param object  Seurat object
#' @param Vsegment  Name of V segment assignment meta.data column
#' @param Jsegment  Name of J segment assignment meta.data column
#' @param Dsegment  Name of D segment assignment meta.data column (not currently implemented)
#' @param clone Vector of clone assignment
#' @param cells Which cells should be included
#'
#' @export
#' @importFrom Seurat FetchData

plot_vdj_circos_seurat <- function(object, Vsegment, Jsegment, Dsegment=NULL, clone, cells=c(), ...){

  data <- Seurat::FetchData(object, vars=c(Vsegment, Jsegment, clone))
  colnames(data) <- c("V","J","clone")

  if(length(cells)>0) data <- data[cells, ]

  return(plot_vdj_circos(Vsegment=data$V, Jsegment=data$J, clone=data$clone, ...))
}

#' Make circos plot into grob
#'
#' A bit of a hack until I find a better way. This allows circos plots to be included in patchwork/plot_grid panels.
#'
#' @param plot_FUN The full circos plot function call
#'
#' @export
#' @importFrom cowplot as_grob

circos_grob <- function(plot_FUN){
  plot <- function() {plot_FUN}

  cowplot::as_grob(plot)
}
