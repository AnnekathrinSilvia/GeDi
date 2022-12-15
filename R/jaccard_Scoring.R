#' Calculate the Jaccard distance
#'
#' Calculate the Jaccard distance for two given genesets.
#'
#' @param a,b A vector of gene names whose interactions should
#'            be scored.
#'
#' @return The Jaccard distance of the two genesets
#' @import dplyr
#' @export
#'
#' @examples
#' a <- c("PDHB", "VARS2")
#' b <- c("IARS2", "PDHA1")
#' c <- calculateJaccard(a, b)
calculateJaccard <- function(a, b){
  if(length(a) == 0 || length(b) == 0){
    return(1)
  }

  set_int <- length(intersect(a, b))
  jaccard <- set_int / (length(a) + length(b) - set_int)

  return(1 - jaccard)
}

#' Get Matrix of Jaccard distances
#'
#' Calculate the Jaccard distance of all combinations of genesets in a given data
#' set of genesets.
#'
#' @param genesets A `list` of genesets (each geneset is represented by a `list`
#'                 of the corresponding genes).
#' @param progress An optional [shiny::Progress()] object to track the progress
#'                 progress of the function in the app.
#'
#' @return A [Matrix::Matrix()] with the pairwise Jaccard distance of each
#'         geneset pair.
#' @export
#' @import parallel
#'
#' @examples
#' genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
#' m <- getJaccardMatrix(genesets)
getJaccardMatrix <- function(genesets, progress = NULL){
  l <- length(genesets)
  if(l == 0){
    return(NULL)
  }
  j <- Matrix::Matrix(0, l, l)
  results <- list()

  n_cores <- parallel::detectCores()
  n_cores <- max(round(n_cores / 2), 1)

  for (k in 1:(l-1)) {
    a <- genesets[[k]]
    if (!is.null(progress)) {
      progress$inc(1 / (l + 1), detail = paste("Scoring geneset number", k))
    }
    results[[k]] <- parallel::mclapply((k+1):l, function(i) {
      b <- genesets[[i]]
      calculateJaccard(a, b)
    }, mc.cores = n_cores)
    j[k,(k+1):l] <- j[(k+1):l, k] <- unlist(results[[k]])
  }

  return(round(j, 2))
}
