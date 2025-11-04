#' Get Matrix of Kappa distances
#'
#' Calculate the Kappa distance of all combinations of genesets in a given data
#' set of genesets. The Kappa distance is normalized to the (0, 1) interval.
#'
#' @param genesets a `list`, A `list` of genesets where each genesets is represented
#'                 by `list` of genes.
#'
#' @return A [Matrix::Matrix()] with Kappa distance rounded to 2 decimal
#'         places.
#' @export
#' @importFrom Matrix Matrix
#' @importFrom proxyC simil
#'
#' @examples
#' #' ## Mock example showing how the data should look like
#' genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
#' m <- getKappaMatrix(genesets)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' genes <- GeDi::prepareGenesetData(macrophage_topGO_example_small)
#' kappa <- getKappaMatrix(genes)
getKappaMatrix <- function(genesets) {
  stopifnot("You must provide at least 2 genesets" = length(genesets) >= 2)
  
  all <- unique(unlist(genesets))
  genesets <- lapply(genesets, function(x) as.numeric(factor(x, levels = all)))
  n <- length(genesets)

  mg <- matrix(0, ncol = length(all), nrow = n)
  for(i in seq_len(n)) {
    mg[i, genesets[[i]]] <- 1
  }
  mg <- as(mg, "sparseMatrix")

  tab <- ncol(mg)
  oab <- proxyC::simil(mg, method = "simple matching")
  m1 <- rowSums(mg)
  m2 <- abs(rowSums(mg - 1))
  aab <- (outer(m1, m1) + outer(m2, m2))/tab/tab
  k <- (oab - aab)/(1 - aab)
  
  k <- round((1-k)/2, 2)
  k@x[!is.finite(k@x)] <- 0

  return(k)
}

