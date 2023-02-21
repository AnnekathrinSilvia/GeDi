#' Get Matrix of Meet-Min distances
#'
#' Calculate the Meet-Min distance of all combinations of genesets in a given
#' data set of genesets.
#'
#' @param genesets A `list` of genesets (each geneset is represented by a `list`
#'                 of the corresponding genes).
#' @param progress An optional [shiny::Progress()] object to track the progress
#'                 of the function in the app.
#' @param n_cores Numerical value indicating the number of cores to be used.
#'                If no value is given, half of the available cores will be
#'                used.
#'
#' @return A [Matrix::Matrix()] with the pairwise Meet-Min distance of each
#'         geneset pair. The matrix is symmetrical with values between 0 and 1,
#'         where 0 indicates the smallest distance (identical genesets) and
#'         1 indicates two disjoint sets.
#' @export
#' @import parallel
#'
#' @examples
#' genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
#' m <- getMeetMinMatrix(genesets, n_cores = 1)
getMeetMinMatrix <- function(genesets, progress = NULL, n_cores = NULL) {
  l <- length(genesets)
  if (l == 0) {
    return(NULL)
  }
  m <- Matrix::Matrix(0, l, l)
  results <- list()

  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores()
    n_cores <- max(round(n_cores / 2), 1)
  }

  for (j in 1:(l - 1)) {
    a <- genesets[[j]]
    if (!is.null(progress)) {
      progress$inc(1 / l, detail = paste("Scoring geneset number", j))
    }
    results[[j]] <- parallel::mclapply((j + 1):l, function(i) {
      b <- genesets[[i]]
      if (length(a) == 0 || length(b) == 0) {
        return(1)
      } else {
        int <- length(intersect(a, b))
        return(1 - (int / min(length(a), length(b))))
      }
    }, mc.cores = n_cores)
    m[j, (j + 1):l] <- m[(j + 1):l, j] <- unlist(results[[j]])
  }
  return(m)
}
