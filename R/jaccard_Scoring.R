#' Calculate the Jaccard distance
#'
#' Calculate the Jaccard distance between two genesets.
#'
#' @param a,b character vector, set of gene identifiers.
#'
#' @return The Jaccard distance of the sets.
#' @export
#'
#' @examples
#' ## Mock example showing how the data should look like
#' a <- c("PDHB", "VARS2")
#' b <- c("IARS2", "PDHA1")
#' c <- calculateJaccard(a, b)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' genes <- GeDi::prepareGenesetData(macrophage_topGO_example_small)
#' jaccard <- calculateJaccard(genes[1], genes[2])
calculateJaccard <- function(a, b) {
  # Calculate the lengths of the input sets
  len_a <- length(a)
  len_b <- length(b)
  # If either set is empty, return a Jaccard distance of 1
  if (len_a == 0 || len_b == 0) {
    return(1)
  }

  # Calculate the size of the intersection between the sets
  set_int <- length(intersect(a, b))
  # Calculate the Jaccard similarity coefficient
  jaccard <- set_int / (len_a + len_b - set_int)
  # Calculate the Jaccard distance by subtracting the Jaccard coefficient from 1
  return(1 - jaccard)
}

#' Get Matrix of Jaccard distances
#'
#' Calculate the Jaccard distance of all combinations of genesets in a given 
#' data set of genesets.
#'
#' @param genesets a `list`, A `list` of genesets where each genesets is 
#'                 represented by `list` of genes.
#' @param progress a [shiny::Progress()] object, Optional progress bar object
#'                 to track the progress of the function (e.g. in a Shiny app).
#' @param BPPARAM A BiocParallel `bpparam` object specifying how parallelization
#'                should be handled. Defaults to [BiocParallel::SerialParam()]
#'
#' @return A [Matrix::Matrix()] with Jaccard distance rounded to 2 decimal
#'         places.
#' @export
#' @importFrom parallel mclapply
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom Matrix Matrix
#'
#' @examples
#' ## Mock example showing how the data should look like
#' genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
#' m <- getJaccardMatrix(genesets)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' genes <- GeDi::prepareGenesetData(macrophage_topGO_example_small)
#' jaccard <-getJaccardMatrix(genes)
getJaccardMatrix <- function(genesets,
                             progress = NULL,
                             BPPARAM = BiocParallel::SerialParam()) {
  # Get the number of gene sets
  l <- length(genesets)
  # If there are no gene sets, return NULL
  if (l == 0) {
    return(NULL)
  }

  # Initialize an empty matrix for storing Jaccard distances
  j <- Matrix(0, l, l)

  # Initialize a list for storing intermediate results
  results <- list()
  # Calculate the Jaccard distance for each pair of gene sets
  for (k in seq_len((l - 1))) {
    a <- genesets[[k]]
    # Update the progress bar if provided
    if (!is.null(progress)) {
      progress$inc(1 / (l + 1), detail = paste("Scoring geneset number", k))
    }
    # Parallelly calculate Jaccard distances for pairs
    results[[k]] <- bplapply((k + 1):l, function(i){
      b <- genesets[[i]]
      calculateJaccard(a, b)
    }, BPPARAM = BPPARAM)
    # Fill the upper and lower triangular sections of the matrix with results
    j[k, (k + 1):l] <- j[(k + 1):l, k] <- unlist(results[[k]])
  }
  # Return the Jaccard distance matrix rounded to 2 decimal places
  return(round(j, 2))
}
