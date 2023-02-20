#' TODO: Title
#'
#' @param genes A `list` of `list` of genes which belong to the genesets in
#'              geneset_names
#' @param geneset_names A vector of geneset names
#' @param start numeric, at which number to start the x-axis, defaults to NULL
#'              meaning the x-axis starts at 0
#' @param end numeric, at which number to end the x-axis, defaults to NULL,
#'        meaning the size of the largest geneset is the end
#' @param binwidth numeric, size of the individual bins, defaults to 6
#' @param color color to use for the bars of the histogram, defaults to #0092AC
#'
#' @return A `ggplot2` histogram with the distribution of the size of the
#'         genesets, i.e. number of genes in the genesets
#'
#' @export
#' @import ggplot2
#'
#' @examples
#' geneset_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
#' genes <- list(c("PDHB", "VARS2"), c("IARS2", "PDHA1"),
#'  c("AAAS", "ABCE1"), c("ABI1", "AAR2"), c("AATF", "AMFR"),
#'  c("BMS1", "DAP3"), c("AURKAIP1", "CHCHD1"), c("IARS2"),
#'  c("AHI1", "ALMS1"))
#'
#'  p <- gs_histogram(genes, geneset_names)
gs_histogram <- function(genes,
                         geneset_names,
                         start = NULL,
                         end = NULL,
                         binwidth = 5,
                         color = "#0092AC"){
  n_genes <- .buildHistogramData(genes, geneset_names, start, end)

  p <- ggplot(n_genes, aes(x = Size)) +
    geom_histogram(binwidth = binwidth, fill = color) +
    theme_bw()

  return(p)
}


#' TODO: Title
#'
#' @param genes A `list` of `list` of genes which belong to the genesets in
#'              geneset_names
#' @param geneset_names A vector of geneset names
#' @param start numeric, at which number to start the x-axis, defaults to NULL
#'              meaning the x-axis starts at 0
#' @param end numeric, at which number to end the x-axis, defaults to NULL,
#'        meaning the size of the largest geneset is the end
#'
#' @return A `data.frame` which maps each geneset to its size (i.e. number of
#'         genes in the genesets), if start and/ or end are given, only those
#'         genesets whose size are larger/smaller than start/end are returned
#'
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
