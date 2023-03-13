#' Calculate interaction score for two genesets
#'
#' The function calculates the interaction score for two given genesets.
#'
#' @param a,b character vector, set of gene identifiers.
#' @param ppi `data.frame`, a `data.frame` object which contains the
#'            Protein-Protein interactions information. The object is expected
#'            to have three columns, `Gene1` and `Gene2` which specify the gene
#'            names of the interacting proteins in no particular order
#'            (symmetric interaction) and a column `combined_score` which is a
#'            numerical value of the strength of the interaction.
#' @param maxInteract numeric, largest value of `ppi`.
#'
#' @return The interaction score of the two given genesets.
#' @export
#' @import dplyr
#'
#' @examples
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
getInteractionScore <- function(a, b, ppi, maxInteract) {
  # get size of set a and b
  len_a <- length(a)
  len_b <- length(b)

  # get the set of genes exclusive to a, b and the intersection of both
  onlya <- setdiff(a, b)
  onlyb <- setdiff(b, a)
  int <- intersect(a, b)

  w <- min(len_a, len_b) / (len_a + len_b)

  # get the size of the sets
  len_int <- length(int)
  len_onlyB <- length(onlyb)
  len_onlyA <- length(onlya)
  lengths <- c(len_a, len_b, len_int, len_onlyA, len_onlyB)

  # check if any set is empty, if so score will be 0
  if (any(lengths == 0)) {
    return(0)
  }

  # calculate sums of interacting genes
  # sumInt <- sum(ppi[which(ppi$Gene1 %in% onlya & ppi$Gene2 %in% int), ]$combined_score)
  # sumOnlyB <- sum(ppi[which(ppi$Gene1 %in% onlya & ppi$Gene2 %in% onlyb), ]$combined_score)

  sumInt <- sum(ppi[ppi$Gene1 %in% onlya & ppi$Gene2 %in% int, "combined_score"])
  sumOnlyB <- sum(ppi[ppi$Gene1 %in% onlya & ppi$Gene2 %in% onlyb, "combined_score"])

  nom <- (w * sumInt) + sumOnlyB
  denom <- maxInteract * (w * len_int + len_onlyB)

  score <- nom / denom

  return(score)
}

#' Calculate local pMM distance
#'
#' Calculate the local pMM distance of two genesets.
#'
#' @param a,b character vector, set of gene identifiers.
#' @param ppi `data.frame`, a `data.frame` object which contains the
#'            Protein-Protein interactions information. The object is expected
#'            to have three columns, `Gene1` and `Gene2` which specify the gene
#'            names of the interacting proteins in no particular order
#'            (symmetric interaction) and a column `combined_score` which is a
#'            numerical value of the strength of the interaction.
#' @param maxInteract numeric, largest value of `ppi`.
#' @param alpha numeric, scaling factor in (0, 1) which indicates how much the
#'              Protein-Protein interactions are weighted into the final score.
#'              Defaults to 1.
#'
#' @return The pMM distance for two sets.
#' @export
#'
#' @examples
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
#' pMM_score <- pMMlocal(a, b, ppi, maxInteract)
pMMlocal <- function(a, b, ppi, maxInteract, alpha = 1) {
  z <- min(length(a), length(b))
  if (z == 0) {
    return(1)
  }

  factor1 <- (length(intersect(a, b))) / z
  interaction_score <- min(
    getInteractionScore(a, b, ppi, maxInteract),
    getInteractionScore(b, a, ppi, maxInteract)
  )
  factor2 <- (alpha / z) * interaction_score

  return(min(factor1 + factor2, 1))
}


#' Calculate the pMM distance
#'
#' Calculate the pMM distance of all combinations of genesets in a given data
#' set of genesets.
#'
#' @param genes `list`, a `list` of genesets (each geneset is represented by
#'               a `list` of the corresponding genes).
#' @param ppi `data.frame`, a `data.frame` object which contains the
#'            Protein-Protein interactions information. The object is expected
#'            to have three columns, `Gene1` and `Gene2` which specify the gene
#'            names of the interacting proteins in no particular order
#'            (symmetric interaction) and a column `combined_score` which is a
#'            numerical value of the strength of the interaction.
#' @param alpha numeric, scaling factor in (0, 1) which indicates how much the
#'              Protein-Protein interactions are weighted into the final score.
#'              Defaults to 1.
#' @param progress [shiny::Progress()] object, optional. To track the progress
#'                 of the function (e.g. in a Shiny app)
#' @param n_cores numeric, number of cores to use for the function.
#'                Defaults to `NULL` in which case the function takes half of
#'                the available cores (see function .detectNumberCores(n_cores)).
#'
#' @return A [Matrix::Matrix()] with the pairwise pMM distance of each
#'         geneset pair.
#' @export
#' @import parallel
#' @import Matrix
#'
#' @examples
#' genes <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
#'
#' ppi <- data.frame(
#'   Gene1 = c("PDHB", "VARS2"),
#'   Gene2 = c("IARS2", "PDHA1"),
#'   combined_score = c(0.5, 0.2)
#' )
#'
#' pMM <- getpMMMatrix(genes, ppi, n_cores = 1)
getpMMMatrix <- function(genes, ppi, alpha = 1, progress = NULL, n_cores = NULL) {
  l <- length(genes)
  if (l == 0) {
    return(NULL)
  }
  stopifnot(!is.null(ppi))
  scores <- Matrix::Matrix(0, l, l)
  maxInteract <- max(ppi$combined_score)

  # determine number of cores to use
  n_cores <- .getNumberCores(n_cores)

  results <- list()

  # calculate pMM distance for each pair of genesets
  for (j in 1:(l - 1)) {
    a <- genes[[j]]
    if (!is.null(progress)) {
      progress$inc(1 / l, detail = paste("Scoring geneset number", j))
    }
    results[[j]] <- parallel::mclapply((j + 1):l, function(i) {
      b <- genes[[i]]
      pMMlocal(a, b, ppi, maxInteract)
    }, mc.cores = n_cores)
    scores[j, (j + 1):l] <- scores[(j + 1):l, j] <- unlist(results[[j]])
  }
  return(scores)

}
