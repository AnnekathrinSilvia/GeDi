#' Create a histogram plot for gene set sizes
#'
#' Create a histogram plot to plot geneset names / identifiers againts their size.
#'
#' @param genesets a `list`, A `list` of genesets where each genesets is represented
#'                 by `list` of genes.
#' @param gs_names character vector, Name / identifier of the genesets in
#'                 `genesets`
#' @param start numeric, Optional, describes the minimum gene set size to
#'              include. Defaults to 0.
#' @param end numeric, Optional, describes the maximum gene set size to include.
#'            Defaults to 0.
#' @param binwidth numeric, Width of histogram bins. Defaults to 5.
#' @param color character, Fill color for histogram bars. Defaults to #0092AC.
#'
#' @return A [ggplot2::ggplot()] plot object.
#'
#' @export
#' @import ggplot2
#'
#' @examples
#' gs_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
#' genesets <- list(
#'   c("PDHB", "VARS2"), c("IARS2", "PDHA1"),
#'   c("AAAS", "ABCE1"), c("ABI1", "AAR2"), c("AATF", "AMFR"),
#'   c("BMS1", "DAP3"), c("AURKAIP1", "CHCHD1"), c("IARS2"),
#'   c("AHI1", "ALMS1")
#' )
#'
#' p <- gs_histogram(genesets, gs_names)
gs_histogram <- function(genesets,
                         gs_names,
                         start = NULL,
                         end = NULL,
                         binwidth = 5,
                         color = "#0092AC") {
  # Check if there are gene sets provided
  stopifnot(length(genesets) > 0)

  # Build a data frame containing gene set sizes
  n_genes <- .buildHistogramData(genesets, gs_names, start, end)

  # Create a histogram plot using ggplot2
  p <- ggplot(n_genes, aes(x = Size)) +
    geom_histogram(binwidth = binwidth, fill = color) +
    theme_bw()

  # Return the histogram plot object
  return(p)
}


#' Prepare data for \code{gs_histogram()}.
#'
#' Prepare the data for the \code{gs_histogram()} by generating a `data.frame`
#' which maps geneset names / identifiers to the size of their size.
#'
#' @param genesets a `list`, A `list` of genesets where each genesets is represented
#'              by `list` of genes.
#' @param gs_names character vector, Name / identifier of the genesets in
#'                 `genesets`
#' @param start numeric, Optional, describes the minimum gene set size to
#'              include. Defaults to 0.
#' @param end numeric, Optional, describes the maximum gene set size to include.
#'            Defaults to 0.
#'
#' @return A `data.frame` mapping geneset names to sizes
#'
.buildHistogramData <- function(genesets,
                                gs_names,
                                start = 0,
                                end = 0) {
  # Get the size of each gene set
  n_genes <- sapply(genesets, length)

  # Create a data frame to store geneset sizes
  n_genes <- as.data.frame(n_genes)
  colnames(n_genes) <- "Size"

  # Check if the number of geneset names matches the number of genesets
  stopifnot(length(gs_names) == length(genesets))

  # Add geneset names to the data frame
  n_genes$Geneset <- gs_names

  # Filter genesets based on the specified start and end sizes (if provided)
  if (!is.null(start) && !is.null(end)) {
    n_genes <- n_genes[(n_genes$Size >= start & n_genes$Size <= end), ]
  }

  # Reorder columns to have geneset names first
  n_genes <- n_genes[, c(2, 1)]

  # Return the final data frame
  return(n_genes)
}
