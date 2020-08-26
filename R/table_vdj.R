table_vdj <- function(vdj, topN=10){
  vdjTable <- vdj %>%
    group_by(chain, v_gene, d_gene, j_gene, c_gene, cdr3) %>%
    summarize(cells=n(),
              reads=sum(reads),
              umis=sum(umis)) %>%
    arrange(-cells, -umis, -reads) %>%
    ungroup() %>%
    mutate(pct=round(cells/sum(cells)*100,2)) %>%
    top_n(topN)

  return(vdjTable)
}
