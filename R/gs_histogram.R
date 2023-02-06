#' Title
#'
#' @param genes
#' @param start
#' @param end
#'
#' @return
#' @export
#' @import ggplot2
#'
#' @examples
gs_histogram <- function(genes,
                         gs_names,
                         start = NULL,
                         end = NULL,
                         binwidth = 5,
                         color = "#0092AC"){
  # get the number of genes per geneset
  n_genes <- .buildHistogramData(genes, gs_names, start, end)

  p <- ggplot(n_genes, aes(x = Size)) +
    geom_histogram(binwidth = binwidth, fill = color) +
    theme_bw()

  return(p)
}


#' Title
#'
#' @param genes
#' @param gs_names
#' @param start
#' @param end
#'
#' @return
#' @export
#'
#' @examples
.buildHistogramData <- function(genes,
                                gs_names,
                                start = NULL,
                                end = NULL){
  n_genes <- sapply(genes, length)
  n_genes <- as.data.frame(n_genes)
  colnames(n_genes) <- "Size"
  n_genes$Geneset <- gs_names

  # filter only for those between start and end
  if(!is.null(start) && !is.null(end)){
    n_genes <- subset(n_genes, Size >= start & Size <= end, select = c(Geneset, Size))
  }

  return(n_genes)
}
