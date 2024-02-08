#' Create a histogram plot for gene set sizes
#'
#' Create a histogram plot to plot geneset names / identifiers against their size.
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
#' ## Mock example showing how the data should look like
#' gs_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
#' genesets <- list(
#'   c("PDHB", "VARS2"), c("IARS2", "PDHA1"),
#'   c("AAAS", "ABCE1"), c("ABI1", "AAR2"), c("AATF", "AMFR"),
#'   c("BMS1", "DAP3"), c("AURKAIP1", "CHCHD1"), c("IARS2"),
#'   c("AHI1", "ALMS1")
#' )
#'
#' p <- gsHistogram(genesets, gs_names)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' genes <- GeDi::getGenes(macrophage_topGO_example_small)
#' p <- gsHistogram(genes, macrophage_topGO_example_small$Genesets)
gsHistogram <- function(genesets,
                         gs_names,
                         start = 0,
                         end = 0,
                         binwidth = 5,
                         color = "#0092AC") {
  # Check if there are gene sets provided
  stopifnot(length(genesets) > 0)

  # Build a data frame containing gene set sizes
  n_genes <- buildHistogramData(genesets, gs_names, start, end)

  # Create a histogram plot using ggplot2
  p <- ggplot(n_genes, aes(x = Size)) +
    geom_histogram(binwidth = binwidth, fill = color) +
    theme_bw()

  # Return the histogram plot object
  return(p)
}


#' Prepare data for \code{gsHistogram()}.
#'
#' Prepare the data for the \code{gsHistogram()} by generating a `data.frame`
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
#' @export
#' @return A `data.frame` mapping geneset names to sizes
#'
#' @examples
#' ## Mock example showing how the data should look like
#' gs_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
#' genesets <- list(
#'   c("PDHB", "VARS2"), c("IARS2", "PDHA1"),
#'   c("AAAS", "ABCE1"), c("ABI1", "AAR2"), c("AATF", "AMFR"),
#'   c("BMS1", "DAP3"), c("AURKAIP1", "CHCHD1"), c("IARS2"),
#'   c("AHI1", "ALMS1")
#' )
#'
#' p <- buildHistogramData(genesets, gs_names)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' genes <- GeDi::getGenes(macrophage_topGO_example_small)
#' p <- buildHistogramData(genes, macrophage_topGO_example_small$Genesets)
buildHistogramData <- function(genesets,
                                gs_names,
                                start = 0,
                                end = 0) {
  # Get the size of each gene set
  n_genes <- vapply(genesets, length, numeric(1))

  # Create a data frame to store geneset sizes
  n_genes <- as.data.frame(n_genes)
  colnames(n_genes) <- "Size"

  # Check if the number of geneset names matches the number of genesets
  stopifnot(length(gs_names) == length(genesets))

  # Add geneset names to the data frame
  n_genes$Geneset <- gs_names

  # Filter genesets based on the specified start and end sizes (if provided)
  if(start > 0){
    n_genes <- n_genes[(n_genes$Size >= start), ]
  }

  if(end > 0){
    n_genes <- n_genes[(n_genes$Size <= end), ]
  }


  # Reorder columns to have geneset names first
  n_genes <- n_genes[, c(2, 1)]

  # Return the final data frame
  return(n_genes)
}
