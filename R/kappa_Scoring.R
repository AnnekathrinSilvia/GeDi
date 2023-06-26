#' Calculate the Kappa distance
#'
#' Calculate the Kappa distance for two given genesets.
#'
#' @param a,b character vector, set of gene identifiers.
#' @param all_genes character vector, list of all (unique) genes available in
#'                  the input data.
#'
#' @return The Kappa distance of the two sets.
#' @export
#'
#' @examples
#' a <- c("PDHB", "VARS2")
#' b <- c("IARS2", "PDHA1")
#' all_genes <- c("PDHB", "VARS2", "IARS2", "PDHA1")
#' c <- calculateKappa(a, b, all_genes)
calculateKappa <- function(a, b, all_genes) {
  n_genes <- length(all_genes)
  stopifnot(n_genes > 0)
  if (length(a) == 0 || length(b) == 0) {
    return(1)
  }

  set_int <- length(intersect(a, b))
  set_none <- sum(!all_genes %in% a & !all_genes %in% b)

  only_a <- sum(all_genes %in% a & !all_genes %in% b)
  only_b <- sum(!all_genes %in% a & all_genes %in% b)

  total <- sum(only_a, only_b, set_int, set_none)

  O <- (set_int + set_none) / total
  E <- (set_int + only_a) * (set_int + only_b) + (only_b + set_none) * (only_a + set_none)
  E <- E / total^2

  kappa <- ((O - E) / (1 - E))

  if (is.nan(kappa)) {
    kappa <- 1
  }
  return(kappa)
}

#' Get Matrix of Kappa distances
#'
#' Calculate the Kappa distance of all combinations of genesets in a given data
#' set of genesets. The Kappa distance is normalized to the (0, 1) interval.
#'
#' @param genesets `list`, a `list` of genesets (each geneset is represented by
#'                 a `list` of the corresponding genes).
#' @param progress [shiny::Progress()] object, optional. To track the progress
#'                 of the function (e.g. in a Shiny app)
#' @param n_cores numeric, number of cores to use for the function.
#'                Defaults to `NULL` in which case the function takes half of
#'                the available cores (see function .detectNumberCores(n_cores)).
#'
#' @return A [Matrix::Matrix()] with the pairwise Kappa distance of each
#'         geneset pair. The matrix is symmetrical with values between 0 and 1,
#'         where 0 indicates the smallest distance (identical genesets) and
#'         1 indicates two disjoint sets.
#' @export
#' @import parallel
#'
#' @examples
#' genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
#' m <- getKappaMatrix(genesets, n_cores = 1)
getKappaMatrix <- function(genesets, progress = NULL, n_cores = NULL) {
  l <- length(genesets)
  if (l == 0) {
    return(NULL)
  }

  # set up parameters
  k <- Matrix::Matrix(0, l, l)
  unique_genes <- unique(unlist(genesets))

  # determine number of cores to use
  n_cores <- .getNumberCores(n_cores)

  results <- list()

  # calculate Kappa distance for each pair of genesets
  for (j in 1:(l - 1)) {
    a <- genesets[[j]]
    if (!is.null(progress)) {
      progress$inc(1 / (l + 1), detail = paste("Scoring geneset number", j))
    }
    results[[j]] <- parallel::mclapply((j + 1):l, function(i) {
      b <- genesets[[i]]
      calculateKappa(a, b, unique_genes)
    }, mc.cores = n_cores)
    k[j, (j + 1):l] <- k[(j + 1):l, j] <- unlist(results[[j]])
  }

  # normalize each value to the (0, 1) interval
  min <- min(k)
  max <- max(k)
  if (!is.null(progress)) {
    progress$inc(1 / (l + 1), detail = "Normalizing Kappa Matrix")
  }

  results <- list()
  for (j in 1:(l - 1)) {
    results[[j]] <- parallel::mclapply((j + 1):l, function(i) {
      return(1 - ((k[j, i] - min) / (max - min)))
    }, mc.cores = n_cores)
    k[j, (j + 1):l] <- k[(j + 1):l, j] <- unlist(results[[j]])
  }

  return(round(k, 2))
}
