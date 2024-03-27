#' Calculate the Kappa distance
#'
#' Calculate the Kappa distance between two genesets.
#'
#' @param a,b character vector, set of gene identifiers.
#' @param all_genes character vector, list of all (unique) genes available in
#'                  the input data.
#'
#' @return The Kappa distance of the sets.
#' @export
#'
#' @examples
#' ## Mock example showing how the data should look like
#' a <- c("PDHB", "VARS2")
#' b <- c("IARS2", "PDHA1")
#' all_genes <- c("PDHB", "VARS2", "IARS2", "PDHA1")
#' c <- calculateKappa(a, b, all_genes)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' genes <- GeDi::getGenes(macrophage_topGO_example_small)
#' c <- calculateKappa(genes[1], genes[2], unique(genes))
calculateKappa <- function(a, b, all_genes) {
  # Get the total number of genes
  n_genes <- length(all_genes)

  # Ensure there are genes to work with
  stopifnot(n_genes > 0)

  # If either set is empty, return Kappa distance of 1
  if (length(a) == 0 || length(b) == 0) {
    return(1)
  }

  # Calculate the size of the intersection between the sets
  set_int <- length(intersect(a, b))

  # Calculate the size of genes not present in either set
  set_none <- sum(!all_genes %in% a & !all_genes %in% b)

  # Calculate the sizes of genes only in one of the sets
  only_a <- sum(all_genes %in% a & !all_genes %in% b)
  only_b <- sum(!all_genes %in% a & all_genes %in% b)

  # Calculate the total number of comparisons
  total <- sum(only_a, only_b, set_int, set_none)

  # Calculate observed agreement (O) and expected agreement (E)
  O <- (set_int + set_none) / total
  E <- (set_int + only_a) * (set_int + only_b) + (only_b + set_none) * (only_a + set_none)
  E <- E / total^2

  # Calculate Cohen's Kappa coefficient
  kappa <- ((O - E) / (1 - E))

  # Handle the case of NaN (e.g., when E is 1)
  if (is.nan(kappa)) {
    kappa <- 1
  }

  # Return the calculated Kappa coefficient
  return(kappa)
}

#' Get Matrix of Kappa distances
#'
#' Calculate the Kappa distance of all combinations of genesets in a given data
#' set of genesets. The Kappa distance is normalized to the (0, 1) interval.
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
#' @return A [Matrix::Matrix()] with Kappa distance rounded to 2 decimal
#'         places.
#' @export
#' @importFrom parallel mclapply
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom Matrix Matrix
#'
#' @examples
#' #' ## Mock example showing how the data should look like
#' genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
#' m <- getKappaMatrix(genesets, n_cores = 1)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' genes <- GeDi::getGenes(macrophage_topGO_example_small)
#' kappa <-getKappaMatrix(genes, n_cores = 1)
getKappaMatrix <- function(genesets, progress = NULL, n_cores = NULL) {
  # Get the number of genesets
  l <- length(genesets)

  # If there are no gene sets, return NULL
  if (l == 0) {
    return(NULL)
  }

  # Initialize an empty matrix for storing Kappa distances
  k <- Matrix(0, l, l)

  # Get the unique genes present across all gene sets
  unique_genes <- unique(unlist(genesets))

  # Initialize a list for storing intermediate results
  results <- list()

  if(Sys.info()["sysname"] == "Windows"){
    # Calculate the Jaccard distance for each pair of gene sets
    for (j in seq_len((l - 1))) {
      a <- genesets[[j]]
      # Update the progress bar if provided
      if (!is.null(progress)) {
        progress$inc(1 / (l + 1), detail = paste("Scoring geneset number", k))
      }
      # Parallelly calculate Jaccard distances for pairs
      results[[j]] <- bplapply((j + 1):l, function(i){
        b <- genesets[[i]]
        calculateKappa(a, b, unique_genes)
      }, BPPARAM = SerialParam())
      # Fill the upper and lower triangular sections of the matrix with results
      k[j, (j + 1):l] <- k[(j + 1):l, j] <- unlist(results[[j]])
    }

    # Normalize each value to the (0, 1) interval
    min <- min(k)
    max <- max(k)
    if (!is.null(progress)) {
      progress$inc(1 / (l + 1), detail = "Normalizing Kappa Matrix")
    }

    # Update the matrix with normalized values
    results <- list()
    for (j in seq_len((l - 1))) {
      results[[j]] <- bplapply((j + 1):l, function(i) {
        return(1 - ((k[j, i] - min) / (max - min)))
      }, BPPARAM = SerialParam())
      k[j, (j + 1):l] <- k[(j + 1):l, j] <- unlist(results[[j]])
    }
  } else{
    # Determine the number of CPU cores to use for parallel processing
    n_cores <- .getNumberCores(n_cores)

    # Calculate Kappa distance for each pair of genesets
    for (j in seq_len((l - 1))) {
      a <- genesets[[j]]
      # Update the progress bar if provided
      if (!is.null(progress)) {
        progress$inc(1 / (l + 1), detail = paste("Scoring geneset number", j))
      }
      # Parallelly calculate Kappa distances for pairs
      results[[j]] <- mclapply((j + 1):l, function(i) {
        b <- genesets[[i]]
        calculateKappa(a, b, unique_genes)
      }, mc.cores = n_cores)
      # Fill the upper and lower triangular sections of the matrix with results
      k[j, (j + 1):l] <- k[(j + 1):l, j] <- unlist(results[[j]])
    }

    # Normalize each value to the (0, 1) interval
    min <- min(k)
    max <- max(k)
    if (!is.null(progress)) {
      progress$inc(1 / (l + 1), detail = "Normalizing Kappa Matrix")
    }

    # Update the matrix with normalized values
    results <- list()
    for (j in seq_len((l - 1))) {
      results[[j]] <- parallel::mclapply((j + 1):l, function(i) {
        return(1 - ((k[j, i] - min) / (max - min)))
      }, mc.cores = n_cores)
      k[j, (j + 1):l] <- k[(j + 1):l, j] <- unlist(results[[j]])
    }
  }




  # Return the normalized Kappa distance matrix rounded to 2 decimal places
  return(round(k, 2))
}
