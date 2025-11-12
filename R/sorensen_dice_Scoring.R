#' Get Matrix of Sorensen-Dice distances
#'
#' Calculate the Sorensen-Dice distance of all combinations of genesets in a
#' given data set of genesets.
#'
#' @param genesets a `list`, A `list` of genesets where each genesets is
#'                 represented by `list` of genes.
#'
#' @return A [Matrix::Matrix()] with Sorensen-Dice distance rounded to 2 decimal
#'         places.
#' @export
#' @importFrom Matrix Matrix
#' @importFrom proxyC simil
#'
#' @examples
#' ## Mock example showing how the data should look like
#' genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
#' m <- getSorensenDiceMatrix(genesets)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' genes <- GeDi::prepareGenesetData(macrophage_topGO_example_small)
#' sd_matrix <- getSorensenDiceMatrix(genes)
getSorensenDiceMatrix <- function(genesets) {
  stopifnot("You must provide at least 2 genesets" = length(genesets) >= 2)
  
  genesets <- genesets
  all <- unique(unlist(genesets))
  genesets <- lapply(genesets, function(x) as.numeric(factor(x, levels = all)))
  n <- length(genesets)

  mg <- matrix(0, ncol = length(all), nrow = n)
  for(i in seq_len(n)) {
    mg[i, genesets[[i]]] <- 1
  }
  mg <- as(mg, "sparseMatrix")

  mat <- proxyC::simil(mg, method = "dice")

  diag(mat) <- 1
  rownames(mat) <- colnames(mat) <- names(genesets)
  mat <- 1 - mat

  return(round(mat, 2))
}
