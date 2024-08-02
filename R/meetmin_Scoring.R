#' Get Matrix of Meet-Min distances
#'
#' Calculate the Meet-Min distance of all combinations of genesets in a given
#' data set of genesets.
#'
#' @param genesets a `list`, A `list` of genesets where each genesets is 
#'                 represented by `list` of genes.
#' @param progress a [shiny::Progress()] object, Optional progress bar object
#'                 to track the progress of the function (e.g. in a Shiny app).
#' @param BPPARAM A BiocParallel `bpparam` object specifying how parallelization
#'                should be handled. Defaults to [BiocParallel::SerialParam()]
#'
#' @return A [Matrix::Matrix()] with Meet-Min distance rounded to 2 decimal
#'         places.
#' @export
#' @importFrom parallel mclapply
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom Matrix Matrix
#'
#' @examples
#' ## Mock example showing how the data should look like
#' genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
#' m <- getMeetMinMatrix(genesets)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' genes <- GeDi::getGenes(macrophage_topGO_example_small)
#' mm <- getMeetMinMatrix(genes)
getMeetMinMatrix <- function(genesets,
                             progress = NULL,
                             BPPARAM = BiocParallel::SerialParam()) {
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
  # Calculate Meet-Min distance for each pair of gene sets
  for (j in seq_len((l - 1))) {
    a <- genesets[[j]]
    # Update the progress bar if provided
    if (!is.null(progress)) {
      progress$inc(1 / l, detail = paste("Scoring geneset number", j))
    }
    # Parallelly calculate Meet-Min distances for pairs
    results[[j]] <- bplapply((j + 1):l, function(i){
      b <- genesets[[i]]
      if (length(a) == 0 || length(b) == 0) {
        return(1)
      } else {
        int <- length(intersect(a, b))
        return(1 - (int / min(length(a), length(b))))
      }
    }, BPPARAM = BPPARAM)
    m[j, (j + 1):l] <- m[(j + 1):l, j] <- unlist(results[[j]])
  }
  # Return the Meet-Min distance matrix rounded to 2 decimal places
  return(round(m, 2))
}
