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
  
  # Calculate the size of sets 
  set_a <- length(a)
  set_b <- length(b)
  
  # If either set is empty, return Kappa distance of 1
  if (set_a == 0 || set_b == 0) {
    return(1)
  }

  # Calculate the size of the intersection between the sets
  set_int <- length(intersect(a, b))
  # Calculate the size of the set union
  set_union <- length(union(a, b))
  # Calculate the complement sets of a and b
  set_compl_a <- length(setdiff(all_genes, a))
  set_compl_b <- length(setdiff(all_genes, b))
  # Calculate the complement set of the intersection of a and b
  set_int_compl <- length(setdiff(all_genes, union(a,b))) 
  
  # Calculate observed agreement (O) and expected agreement (E)
  O <- (set_int + set_int_compl) / n_genes
  E <- (set_a * set_b + set_compl_a * set_compl_b)
  E <- E / (n_genes^2)

  # Calculate Cohen's Kappa coefficient
  kappa <- ((O - E) / (1 - E))
  
  # Handle the case of NaN (i.e E = 1)
  if(is.nan(kappa)){
    # This can happen if both sets are equal. Then kappa should be 0
    if(setequal(a, b)){
      return(0)
    }else{
      # This can also happen when both sets are completely disjoint, then kappa
      # should be -1 and resulting in a distance of 2
      kappa <- -1
    }
  }
  # Calculate the distance as 1 - kappa
  kappa <- 1 - kappa

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
#' @param BPPARAM A BiocParallel `bpparam` object specifying how parallelization
#'                should be handled. Defaults to [BiocParallel::SerialParam()]
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
#' m <- getKappaMatrix(genesets)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' genes <- GeDi::getGenes(macrophage_topGO_example_small)
#' kappa <- getKappaMatrix(genes)
getKappaMatrix <- function(genesets,
                           progress = NULL,
                           BPPARAM = BiocParallel::SerialParam()) {
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
  # Calculate the Jaccard distance for each pair of gene sets
  for (j in seq_len((l - 1))) {
    a <- genesets[[j]]
    # Update the progress bar if provided
    if (!is.null(progress)) {
      progress$inc(1 / (l + 1), detail = paste("Scoring geneset number", j))
    }
    # Parallelly calculate Jaccard distances for pairs
    results[[j]] <- bplapply((j + 1):l, function(i){
      b <- genesets[[i]]
      calculateKappa(a, b, unique_genes)
    }, BPPARAM = BPPARAM)
    # Fill the upper and lower triangular sections of the matrix with results
    k[j, (j + 1):l] <- k[(j + 1):l, j] <- unlist(results[[j]])
  }

  # Normalize each value to the (0, 1) interval by dividing the matrix by 2
  k <- k / 2
 
  # Return the normalized Kappa distance matrix rounded to 2 decimal places
  return(round(k, 2))
}
