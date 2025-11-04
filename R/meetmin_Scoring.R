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
#' genes <- GeDi::prepareGenesetData(macrophage_topGO_example_small)
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

### IMPLEM: this one uses functions from fastmatch ------
getMeetMinMatrix_fm <- function(genesets,
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
        # int <- length(intersect(a, b))

        int <- sum(a %fin% b)

        return(1 - (int / min(length(a), length(b))))
      }
    }, BPPARAM = BPPARAM)
    m[j, (j + 1):l] <- m[(j + 1):l, j] <- unlist(results[[j]])
  }
  # Return the Meet-Min distance matrix rounded to 2 decimal places
  return(round(m, 2))
}


### IMPLEM: this one is possibly just a more efficient rework of the original  ------
meetmin_matrix <- function(genesets) {
  n <- length(genesets)

  M <- matrix(0, n, n, dimnames = list(names(genesets), names(genesets)))

  for (i in seq_len(n)) {
    a <- unique(genesets[[i]])
    for (j in seq(i, n)) {
      b <- unique(genesets[[j]])
      inter <- sum(a %in% b)
      denom <- min(length(a), length(b))
      sim <- if (denom > 0) inter / denom else 0
      M[i, j] <- 1 - sim
      M[j, i] <- M[i, j]
    }
  }
  round(M,2)
}


### IMPLEM: benchmarking the meet min implementations ------
## bench::mark(
##   getMeetMinMatrix(genes),
##   getMeetMinMatrix_fm(genes),
##   meetmin_matrix(genes),
##   check = FALSE,
##   iterations = 100
## )
##
## all.equal(getMeetMinMatrix(genes), meetmin_matrix(genes))
## summary(as.vector(
##   as.matrix(getMeetMinMatrix(genes)) - as.matrix(meetmin_matrix(genes)))
## )
##
## pheatmap::pheatmap(getMeetMinMatrix(genes))
## pheatmap::pheatmap(getMeetMinMatrix_fm(genes))
## pheatmap::pheatmap(meetmin_matrix(genes))
