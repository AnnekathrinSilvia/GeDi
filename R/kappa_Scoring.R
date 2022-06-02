#' Calculate the Kappa distance for two genesets
#'
#' @param a,b individual genesets as two vectors of genes
#' @param all_genes the list of all unique genes in all input genesets
#'
#' @return the Kappa distance of geneset a and b
#' @import dplyr
#' @export
#'
#' @examples
#' a <- c("LARS", "BCAT1")
#' b <- c("BCAT1", "IARS")
#' all_genes <- c("BCAT1", "IARS",  "LARS2")
#' c <- calculateKappa(a, b, all_genes)
calculateKappa <- function(a, b, all_genes){
  set_int <- length(intersect(a, b))
  l <- length(all_genes)

  only_a <- sum(all_genes %in% a & !all_genes %in% b)
  only_b <- sum(!all_genes %in% a & all_genes %in% b)

  background <- l - sum(only_a, only_b, set_int)

  O <- (set_int + background) / l
  E <- (set_int + only_a) * (set_int + only_b) + (only_b + background) * (only_a + background)
  E <- E / l^2

  kappa <- ((O-E) / (1-E))

  if(is.nan(kappa)){
    kappa <- 0
  }
  return(abs(kappa))
}

#' Calculate the Kappa distance of all geneset combinations
#'
#' @param genesets a list of the genesets to score (where the genesets are vectors of genes)
#' @param geneset_names a list of the names/ids of the genesets
#'
#' @return a [Matrix::Matrix()] with the pairwise Kappa distance of each geneset pair
#' @export
#'
#' @examples
#' genesets <- list(c("PDHB", "VARS2", "IARS2"), c("IARS2", "PDHA2"))
#' geneset_names <- list("A", "B")
#' m <- getKappaMatrix(genesets, geneset_names)
getKappaMatrix <- function(genesets, geneset_names){
  l <- length(genesets)
  if(l == 0){
    return(-1)
  }
  k <- Matrix::Matrix(1, l, l)
  unique_genes <- unique(unlist(genesets))

  for(i in 1:(l - 1)){
    a <- unlist(genesets[i])
    for(j in (i+1):l){
      b <- unlist(genesets[j])
      k[i, j] <- k[j, i] <- calculateKappa(a, b, unique_genes)
    }
  }
  rownames(k) <- geneset_names
  colnames(k) <- geneset_names
  return(k)
}
