#' Calculate the Sorensen-Dice distance
#'
#' Calculate the Sorensen-Dice distance between two genesets.
#'
#' @param a,b character vector, set of gene identifiers.
#'
#' @return The Sorensen-Dice distance of the sets.
#' @export
#'
#' @examples
#' a <- c("PDHB", "VARS2")
#' b <- c("IARS2", "PDHA1")
#' c <- calculateSorensenDice(a, b)
calculateSorensenDice <- function(a, b) {
  # Calculate the lengths of the input sets
  len_a <- length(a)
  len_b <- length(b)

  # If either set is empty, return a Sorensen-Dice distance of 1
  if (len_a == 0 || len_b == 0) {
    return(1)
  }

  # Calculate the size of the intersection between the sets
  set_int <- length(intersect(a, b))

  # Calculate the Sorensen-Dice similarity coefficient
  sorensenDice <- 2 * set_int / (len_a + len_b)

  # Calculate the Sorensen-Dice distance by subtracting the coefficient from 1
  return(1 - sorensenDice)
}

#' Get Matrix of Sorensen-Dice distances
#'
#' Calculate the Sorensen-Dice distance of all combinations of genesets in a given data
#' set of genesets.
#'
#' @param genesets a `list`, A `list` of genesets where each genesets is represented
#'                 by `list` of genes.
#' @param progress a [shiny::Progress()] object, Optional progress bar object
#'                 to track the progress of the function (e.g. in a Shiny app).
#' @param n_cores numeric, Optional number of CPU cores to use for parallel
#'                processing. Defaults to `NULL` in which case the function
#'                takes half of the available cores (see
#'                \code{.detectNumberCores()} function).
#'
#' @return A [Matrix::Matrix()] with Sorensen-Dice distance rounded to 2 decimal
#'         places.
#' @export
#' @importFrom parallel mclapply
#' @importFrom Matrix Matrix
#'
#' @examples
#' genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
#' m <- getSorensenDiceMatrix(genesets, n_cores = 1)
getSorensenDiceMatrix <- function(genesets, progress = NULL, n_cores = NULL) {
  # Get the number of gene sets
  l <- length(genesets)

  # If there are no gene sets, return NULL
  if (l == 0) {
    return(NULL)
  }

  # Initialize an empty matrix for storing Sorensen-Dice distances
  s <- Matrix::Matrix(0, l, l)

  # Initialize a list for storing intermediate results
  results <- list()

  # Determine the number of CPU cores to use for parallel processing
  n_cores <- .getNumberCores(n_cores)

  # Calculate the Sorensen-Dice distance for each pair of gene sets
  for (k in seq_len((l - 1))) {
    a <- genesets[[k]]
    # Update the progress bar if provided
    if (!is.null(progress)) {
      progress$inc(1 / (l + 1), detail = paste("Scoring geneset number", k))
    }
    # Parallelly calculate Sorensen-Dice distances for pairs
    results[[k]] <- mclapply((k + 1):l, function(i) {
      b <- genesets[[i]]
      calculateSorensenDice(a, b)
    }, mc.cores = n_cores)
    # Fill the upper and lower triangular sections of the matrix with results
    s[k, (k + 1):l] <- s[(k + 1):l, k] <- unlist(results[[k]])
  }

  # Return the Sorensen-Dice distance matrix rounded to 2 decimal places
  return(round(s, 2))
}
