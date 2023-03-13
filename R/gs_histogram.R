#' Generate a histogram of geneset size
#'
#' Generate a histogram on the size of genesets (i.e. the number of genes
#' associated with the geneset)
#'
#' @param genes `list`, a `list` of `list` of genes which belong to the genesets
#'               in gs_names
#' @param gs_names character vector, names/identifiers of genesets
#' @param start numeric, at which number to start the x-axis. Defaults to NULL
#'              meaning the x-axis starts at 0.
#' @param end numeric, at which number to end the x-axis. Defaults to NULL,
#'        meaning the size of the largest geneset is the end.
#' @param binwidth numeric, size of the individual bins. Defaults to 5.
#' @param color character, color to use for the bars of the histogram.
#'              Defaults to #0092AC.
#'
#' @return A `ggplot2` histogram with the distribution of the size of the
#'         genesets, i.e. number of genes in the genesets.
#'
#' @export
#' @import ggplot2
#'
#' @examples
#' gs_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
#' genes <- list(
#'   c("PDHB", "VARS2"), c("IARS2", "PDHA1"),
#'   c("AAAS", "ABCE1"), c("ABI1", "AAR2"), c("AATF", "AMFR"),
#'   c("BMS1", "DAP3"), c("AURKAIP1", "CHCHD1"), c("IARS2"),
#'   c("AHI1", "ALMS1")
#' )
#'
#' p <- gs_histogram(genes, gs_names)
gs_histogram <- function(genes,
                         gs_names,
                         start = NULL,
                         end = NULL,
                         binwidth = 5,
                         color = "#0092AC") {
  # set up data.frame for histogram
  stopifnot(length(genes) > 0)
  n_genes <- .buildHistogramData(genes, gs_names, start, end)

  # set up histogram plot
  p <- ggplot(n_genes, aes(x = Size)) +
    geom_histogram(binwidth = binwidth, fill = color) +
    theme_bw()

  return(p)
}


#' Generate a data.frame of histogram data
#'
#' Generate a data.frame for a histogram which maps a geneset identifier to the
#' size of the geneset.
#'
#' @param genes `list`, a `list` of `list` of genes which belong to the genesets
#'              in gs_names
#' @param gs_names character vector, names/identifiers of genesets
#' @param start numeric, at which number to start the x-axis. Defaults to NULL
#'              meaning the x-axis starts at 0.
#' @param end numeric, at which number to end the x-axis. Defaults to NULL,
#'        meaning the size of the largest geneset is the end.
#'
#' @return A `data.frame` which maps each geneset to its size (i.e. number of
#'         genes in the genesets). Ff start and/or end are given, only those
#'         genesets whose size is larger/smaller than start/end are returned.
#'
.buildHistogramData <- function(genes,
                                gs_names,
                                start = NULL,
                                end = NULL) {
  # get size of each geneset
  n_genes <- sapply(genes, length)
  # build up data.frame
  n_genes <- as.data.frame(n_genes)
  colnames(n_genes) <- "Size"
  stopifnot(length(gs_names) == length(genes))
  n_genes$Geneset <- gs_names

  # filter only for those between start and end (if given)
  if (!is.null(start) && !is.null(end)) {
    n_genes <- n_genes[(n_genes$Size >= start & n_genes$Size <= end), ]
  }
  n_genes <- n_genes[, c(2, 1)]

  return(n_genes)
}
