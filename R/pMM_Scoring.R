#' Sums up the protein-protein interacions between genes in two genesets
#'
#' @param x,y genesets
#' @param ppi a Protein-Protein interaction matrix
#'
#' @return sum of the Protein-Protein interactions in the two genesets
#'
#' @examples
#' x <- c(1, 2, 3)
#' y <- c(3, 1)
#' ppi <- Matrix::Matrix(0.5, 3, 3)
#' s <- sumInteraction(x, y, ppi)
sumInteraction <- function(x, y, ppi) {
  if (length(x) ==  0 | length(y) == 0) {
    return(0)
  } else{
    return(sum(ppi[x, y]))
  }
}

#' Calculates the interaction score for two genesets
#'
#' @param a,b genesets
#' @param ai,bi indexed genesets
#' @param ppi a Protein-Protein interaction matrix
#' @param maxInteract the maximum value in the Protein-Protein interaction matrix
#'
#' @return the interaction score for two genesets
#'
#' @examples
#' \dontrun{
#' a <- c("PDHB", "VARS2", "IARS2")
#' ai <- c(1, 2, 3)
#' b <- c("IARS2", "PDHA2")
#' bi <- c(4, 3)
#' ppi <- Matrix::Matrix(0.5, 4, 4)
#' maxInteract <- 0.5
#' s <- getInteractionScore(a, ai, b, bi, ppi, maxInteract)}
getInteractionScore <- function(a, ai, b, bi, ppi, maxInteract) {
  if(length(a) == 0 || length(b) == 0 || length(ai) == 0 || length(bi) == 0){
    return(-1)
  }
  onlya <- setdiff(ai, bi)
  onlyb <- setdiff(bi, ai)
  int <- intersect(ai, bi)

  w <- min(length(a), length(b)) / (length(a) + length(b))
  intlength <- length(intersect(a, b))
  onlyblength <- length(setdiff(b, a))

  sumInt <- sumInteraction(onlya, int, ppi)
  sumOnlyb <- sumInteraction(onlya, onlyb, ppi)

  nom <- (w * sumInt) + sumOnlyb
  denom <- maxInteract * (w * intlength + onlyblength)

  return(nom / denom)
}

#' Calculates the pMM distance for two genesets
#'
#' @param a,b genesets
#' @param ai,bi indexed genesets
#' @param alpha scaling factor in (0, 1); indicates how much the Protein-Protein interactions are weighted into the distance
#' @param ppi a Protein-Protein interaction matrix
#' @param maxInteract the maximum value of the Protein-Protein interaction matrix
#'
#' @return pMM distance for two genesets
#' @export
#'
#' @examples
#' a <- c("PDHB", "VARS2", "IARS2")
#' ai <- c(1, 2, 3)
#' b <- c("IARS2", "PDHA2")
#' bi <- c(4, 3)
#' ppi <- Matrix::Matrix(0.5, 4, 4)
#' maxInteract <- 0.5
#' s <- pMMlocal(a, ai, b, bi, ppi, maxInteract)
pMMlocal <- function(a, ai, b, bi, alpha, ppi, maxInteract) {
  z <- min(length(a), length(b))

  factor1 <- (length(intersect(a, b))) / z
  factor2 <- (alpha / z) * getInteractionScore(a, ai, b, bi, ppi, maxInteract)

  return(min(factor1 + factor2, 1))
}


#' Calculate the pMM distance of all geneset combinations
#'
#' @param genes a list of  vectors of genes (the genesets)
#' @param ppi a Protein-Protein interaction matrix
#' @param geneset_names a list of names/ids of the genesets
#' @param alpha scaling factor in (0, 1); indicates how much the Protein-Protein interactions are weighted into the distance
#'
#' @return a [Matrix::Matrix()] with the pairwise pMM distance of each geneset pair
#' @export
#'
#' @examples
#' genes <- list(c("PDHB", "VARS2", "IARS2"), c("IARS2", "PDHA2"))
#' geneset_names <- list("A", "B")
#' ppi <- Matrix::Matrix(0.5, 4, 4)
#' s <- getpMMMatrix(genes, ppi, geneset_names)
getpMMMatrix <- function(genes, ppi, geneset_names, alpha = 1){
  genes_indexed <- getGenesIndexed(genes, ppi)

  l <- length(genes)
  if(l == 0){
    return(-1)
  }
  scores <- Matrix::Matrix(0, l, l)
  maxInteract <- max(ppi)
  for (i in 1:(l - 1)) {
    a <- genes[[i]]
    ai <- genes_indexed[[i]]
    for (j in (i + 1):l) {
      b <- genes[[j]]
      bi <- genes_indexed[[j]]
      pmm <- min(
        pMMlocal(a, ai, b, bi, alpha, ppi, maxInteract),
        pMMlocal(b, bi, a, ai, alpha, ppi, maxInteract)
      )
      scores[i, j] <- scores[j, i] <- pmm
    }
  }

  rownames(scores) <- geneset_names
  colnames(scores) <- geneset_names
  return(scores)
}
