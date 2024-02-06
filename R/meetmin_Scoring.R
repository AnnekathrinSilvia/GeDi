#' Get Matrix of Meet-Min distances
#'
#' Calculate the Meet-Min distance of all combinations of genesets in a given
#' data set of genesets.
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
#' @return A [Matrix::Matrix()] with Meet-Min distance rounded to 2 decimal
#'         places.
#' @export
#' @importFrom parallel mclapply
#' @importFrom Matrix Matrix
#'
#' @examples
#' ## Mock example showing how the data should look like
#' genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
#' m <- getMeetMinMatrix(genesets, n_cores = 1)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' genes <- GeDi::getGenes(macrophage_topGO_example_small)
#' mm <- getMeetMinMatrix(genes, n_cores = 1)
getMeetMinMatrix <- function(genesets, progress = NULL, n_cores = NULL) {
  # Get the number of genesets
  l <- length(genesets)

  # If there are no genesets, return NULL
  if (l == 0) {
    return(NULL)
  }

  # Initialize an empty matrix for storing Meet-Min distances
  m <- Matrix(0, l, l)

  # Initialize a list for storing intermediate results
  results <- list()

  # Determine the number of CPU cores to use for parallel processing
  n_cores <- .getNumberCores(n_cores)

  # Calculate Meet-Min distance for each pair of gene sets
  for (j in seq_len((l - 1))) {
    a <- genesets[[j]]
    # Update the progress bar if provided
    if (!is.null(progress)) {
      progress$inc(1 / l, detail = paste("Scoring geneset number", j))
    }
    # Parallelly calculate Meet-Min distances for pairs
    results[[j]] <- mclapply((j + 1):l, function(i) {
      b <- genesets[[i]]
      if (length(a) == 0 || length(b) == 0) {
        return(1)
      } else {
        int <- length(intersect(a, b))
        return(1 - (int / min(length(a), length(b))))
      }
    }, mc.cores = n_cores)
    # Fill the upper and lower triangular sections of the matrix with results
    m[j, (j + 1):l] <- m[(j + 1):l, j] <- unlist(results[[j]])
  }

  # Return the Meet-Min distance matrix rounded to 2 decimal places
  return(round(m, 2))
}
