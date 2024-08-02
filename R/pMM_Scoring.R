#' Calculate interaction score for two genesets
#'
#' The function calculates an interaction score between two sets of genes based
#' on a protein-protein interaction network.
#'
#' @param a,b character vector, set of gene identifiers.
#' @param ppi a `data.frame`, Protein-protein interaction (PPI) network data 
#'            frame. The object is expected to have three columns, `Gene1` and 
#'            `Gene2` which specify the gene names of the interacting proteins 
#'            in no particular order (symmetric interaction) and a column
#'            `combined_score` which is a numerical value of the strength of
#'            the interaction.
#' @param maxInteract numeric, Maximum interaction value in the PPI.
#'
#' @references
#' See https://doi.org/10.1186/s12864-019-5738-6 for details on the original
#' implementation.
#'
#' @return Interaction score between the two gene sets.
#' @export
#'
#' @examples
#' ## Mock example showing how the data should look like
#' a <- c("PDHB", "VARS2", "IARS2")
#' b <- c("IARS2", "PDHA1")
#'
#' ppi <- data.frame(
#'   Gene1 = c("PDHB", "VARS2", "IARS2"),
#'   Gene2 = c("IARS2", "PDHA1", "CD3"),
#'   combined_score = c(0.5, 0.2, 0.1)
#' )
#' maxInteract <- max(ppi$combined_score)
#'
#' interaction <- getInteractionScore(a, b, ppi, maxInteract)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' genes <- GeDi::getGenes(macrophage_topGO_example_small)
#' data(ppi_macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' maxInteract <- max(ppi_macrophage_topGO_example_small$combined_score)
#'
#' interaction <- getInteractionScore(genes[1], genes[2], ppi, maxInteract)
getInteractionScore <- function(a, b, ppi, maxInteract) {
  # Get the size of sets a and b
  len_a <- length(a)
  len_b <- length(b)
  # Get the sets of genes exclusive to a, b, and the intersection of both
  onlya <- setdiff(a, b)
  onlyb <- setdiff(b, a)
  int <- intersect(a, b)

  # Calculate weight based on the sizes of the sets a and b
  w <- min(len_a, len_b) / (len_a + len_b)
  # Get the size of the intersection and exclusive sets
  len_int <- length(int)
  len_onlyB <- length(onlyb)
  len_onlyA <- length(onlya)
  lengths <- c(len_a, len_b, len_int, len_onlyA, len_onlyB)

  # Check if any set is empty; if so, the score will be 0
  if (any(lengths == 0)) {
    return(0)
  }
  # Calculate sums of interacting genes
  sumInt <-
    sum(ppi[ppi$Gene1 %in% onlya &
              ppi$Gene2 %in% int, "combined_score"])
  sumOnlyB <-
    sum(ppi[ppi$Gene1 %in% onlya &
              ppi$Gene2 %in% onlyb, "combined_score"])
  # Calculate the numerator and denominator for the interaction score
  nom <- (w * sumInt) + sumOnlyB
  denom <- maxInteract * (w * len_int + len_onlyB)
  # Calculate the interaction score
  score <- nom / denom

  # Return the calculated interaction score
  return(score)
}

#' Calculate local pMM distance
#'
#' Calculate the local pMM distance of two genesets.
#'
#' @param a,b character vector, set of gene identifiers.
#' @param ppi a `data.frame`, Protein-protein interaction (PPI) network data 
#'            frame. The object is expected to have three columns, `Gene1` and 
#'            `Gene2` which specify the gene names of the interacting proteins 
#'            in no particular order (symmetric interaction) and a column
#'            `combined_score` which is a numerical value of the strength of
#'            the interaction.
#' @param maxInteract numeric, Maximum interaction value in the PPI.
#' @param alpha numeric, Scaling factor for controlling the influence of the
#'              interaction score. Defaults to 1.
#'
#' @references
#' See https://doi.org/10.1186/s12864-019-5738-6 for details on the original
#' implementation.
#'
#' @return The pMMlocal score between the two gene sets.
#' @export
#'
#' @examples
#' ## Mock example showing how the data should look like
#' a <- c("PDHB", "VARS2")
#' b <- c("IARS2", "PDHA1")
#'
#' ppi <- data.frame(
#'   Gene1 = c("PDHB", "VARS2"),
#'   Gene2 = c("IARS2", "PDHA1"),
#'   combined_score = c(0.5, 0.2)
#' )
#' maxInteract <- max(ppi$combined_score)
#'
#' pMM_score <- pMMlocal(a, b, ppi, alpha = 1,  maxInteract)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' genes <- GeDi::getGenes(macrophage_topGO_example_small)
#' data(ppi_macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' maxInteract <- max(ppi_macrophage_topGO_example_small$combined_score)
#'
#' pMMlocal <- pMMlocal(genes[1], genes[2], ppi, alpha = 1,  maxInteract)
pMMlocal <- function(a, b, ppi, maxInteract, alpha) {
  # Get the minimum size of sets a and b
  z <- min(length(a), length(b))
  # If either set is empty, return pMMlocal score of 1
  if (z == 0) {
    return(1)
  }

  # Calculate factor1 as the proportion of intersection size to z
  factor1 <- (length(intersect(a, b))) / z
  if(alpha == 0){
    return(min(factor1, 1))
  }
  
  # Calculate the interaction score as the minimum of two interaction scores
  interaction_score <- min(
    getInteractionScore(a, b, ppi, maxInteract),
    getInteractionScore(b, a, ppi, maxInteract)
  )
  # Calculate factor2 as a scaled interaction_score using alpha and z
  factor2 <- (alpha / z) * interaction_score
  
  # Return the minimum of the sum of factor1 and factor2, and 1
  return(min(factor1 + factor2, 1))
}


#' Calculate the pMM distance
#'
#' Calculate the pMM distance of all combinations of genesets in a given data
#' set of genesets.
#'
#' @param genesets a `list`, A `list` of genesets where each genesets is 
#'                 represented by `list` of genes.
#' @param ppi a `data.frame`, Protein-protein interaction (PPI) network data 
#'            frame. The object is expected to have three columns, `Gene1` and 
#'            `Gene2` which specify the gene names of the interacting proteins 
#'            in no particular order (symmetric interaction) and a column
#'            `combined_score` which is a numerical value of the strength of
#'            the interaction.
#' @param alpha numeric, Scaling factor for controlling the influence of the
#'              interaction score. Defaults to 1.
#' @param progress a [shiny::Progress()] object, Optional progress bar object
#'                 to track the progress of the function (e.g. in a Shiny app).
#' @param BPPARAM A BiocParallel `bpparam` object specifying how parallelization
#'                should be handled. Defaults to [BiocParallel::SerialParam()]
#'
#' @return A [Matrix::Matrix()] with pMM distance rounded to 2 decimal places.
#'
#' @references
#' See https://doi.org/10.1186/s12864-019-5738-6 for details on the original
#' implementation.
#'
#'
#' @export
#' @importFrom parallel mclapply
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom Matrix Matrix
#'
#' @examples
#' ## Mock example showing how the data should look like
#' genesets <- list(c("PDHB", "VARS2"), c("IARS2", "PDHA1"))
#'
#' ppi <- data.frame(
#'   Gene1 = c("PDHB", "VARS2"),
#'   Gene2 = c("IARS2", "PDHA1"),
#'   combined_score = c(0.5, 0.2)
#' )
#'
#' pMM <- getpMMMatrix(genesets, ppi)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' genes <- GeDi::getGenes(macrophage_topGO_example_small)
#' data(ppi_macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#'
#' pMM <- getpMMMatrix(genes, ppi)
getpMMMatrix <- function(genesets,
                         ppi,
                         alpha = 1,
                         progress = NULL,
                         BPPARAM = BiocParallel::SerialParam()) {
  # Get the number of genesets
  l <- length(genesets)
  # If there are no genesets, return NULL
  if (l == 0) {
    return(NULL)
  }
  # Ensure that the protein-protein interaction (PPI) network is provided
  stopifnot(!is.null(ppi))

  # Initialize an empty matrix for storing pMM distances
  scores <- Matrix(0, l, l)
  # Get the maximum interaction score in the PPI network
  maxInteract <- max(ppi$combined_score)
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
    results[[j]] <- bplapply((j + 1):l, function(i) {
      b <- genesets[[i]]
      pMMlocal(a, b, ppi, maxInteract, alpha)
    }, BPPARAM = BPPARAM)
    scores[j, (j + 1):l] <-
      scores[(j + 1):l, j] <- unlist(results[[j]])
  }
  # Return the pMM distance matrix rounded to 2 decimal places
  return(round(scores, 2))
}
