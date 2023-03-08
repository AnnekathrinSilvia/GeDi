#' Calculate the Jaccard distance
#'
#' Calculate the Jaccard distance for two given genesets.
#'
#' @param a,b character vector, set of gene identifiers.
#'
#' @return The Jaccard distance of the two sets.
#' @import dplyr
#' @export
#'
#' @examples
#' a <- c("PDHB", "VARS2")
#' b <- c("IARS2", "PDHA1")
#' c <- calculateJaccard(a, b)
calculateJaccard <- function(a, b) {
  len_a <- length(a)
  len_b <- length(b)
  if (len_a == 0 || len_b == 0) {
    return(1)
  }

  set_int <- length(intersect(a, b))
  jaccard <- set_int / (len_a + len_b - set_int)

  return(1 - jaccard)
}

#' Get Matrix of Jaccard distances
#'
#' Calculate the Jaccard distance of all combinations of genesets in a given data
#' set of genesets.
#'
#' @param genesets `list`, a `list` of genesets (each geneset is represented
#'                 by a `list` of the corresponding genes).
#' @param progress [shiny::Progress()] object, optional. To track the progress
#'                 of the function (e.g. in a Shiny app).
#' @param n_cores numeric, number of cores to use for the function.
#'                Defaults to `NULL` in which case the function takes half of
#'                the available cores (see function .detectNumberCores(n_cores)).
#'
#' @return A [Matrix::Matrix()] with the pairwise Jaccard distance of each
#'         geneset pair. The matrix is symmetrical with values between 0 and 1,
#'         where 0 indicates the smallest distance (identical genesets) and
#'         1 indicates two disjoint sets.
#' @export
#' @import parallel
#'
#' @examples
#' genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
#' m <- getJaccardMatrix(genesets, n_cores = 1)
getJaccardMatrix <- function(genesets, progress = NULL, n_cores = NULL) {
  l <- length(genesets)
  if (l == 0) {
    return(NULL)
  }
  # set up parameters
  j <- Matrix::Matrix(0, l, l)
  results <- list()

  # determine number of cores to use
  n_cores <- .getNumberCores(n_cores)

  # calculate the Jaccard distance for each pair of genesets
  for (k in 1:(l - 1)) {
    a <- genesets[[k]]
    if (!is.null(progress)) {
      progress$inc(1 / (l + 1), detail = paste("Scoring geneset number", k))
    }
    results[[k]] <- parallel::mclapply((k + 1):l, function(i) {
      b <- genesets[[i]]
      calculateJaccard(a, b)
    }, mc.cores = n_cores)
    j[k, (k + 1):l] <- j[(k + 1):l, k] <- unlist(results[[k]])
  }

  return(round(j, 2))
}
