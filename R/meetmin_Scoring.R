#' Get Matrix of Meet-Min distances
#'
#' Calculate the Meet-Min distance of all combinations of genesets in a given
#' data set of genesets.
#'
#' @param genesets a `list`, A `list` of genesets where each genesets is
#'                 represented by `list` of genes.
#'
#' @return A [Matrix::Matrix()] with Meet-Min distance rounded to 2 decimal
#'         places.
#' @export
#' @importFrom Matrix Matrix
#'
#' @examples
#' ## Mock example showing how the data should look like
#' genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
#' m <- getMeetMinMatrix(genesets)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' genes <- GeDi::prepareGenesetData(macrophage_topGO_example_small)
#' mm <- getMeetMinMatrix(genes)
getMeetMinMatrix <- function(genesets) {
  stopifnot("You must provide at least 2 genesets" = length(genesets) >= 2)
  
  n <- length(genesets)

  M <- matrix(0, n, n, dimnames = list(names(genesets), names(genesets)))

  for (i in seq_len(n)) {
    a <- unique(genesets[[i]])
    for (j in seq(i, n)) {
      b <- unique(genesets[[j]])
      inter <- sum(a %in% b)
      denom <- min(length(a), length(b))
      sim <- if (denom > 0) inter / denom else 0
      M[i, j] <- 1 - sim
      M[j, i] <- M[i, j]
    }
  }
  return(round(M,2))
}

