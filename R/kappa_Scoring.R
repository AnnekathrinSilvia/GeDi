#' Calculate the Kappa distance for two genesets
#'
#' @param a,b individual genesets
#' @param l the number of unique genes in all input genesets
#'
#' @return the Kappa distance of geneset a and b
#' @export
#'
#' @examples
#' a <- c("LARS", "LARS2", "BCAT1")
#' b <- c("BCAT1", "IARS")
#' l <- 4
#' c <- calculateKappa(a, b, l)
calculateKappa <- function(a, b, l){
  set_int <- length(intersect(a, b))

  only_a <- length(setdiff(a, b))
  only_b <- length(setdiff(b, a))

  background <- l - (only_a + only_b + set_int)

  O <- (set_int + background) / l
  E <- (set_int + only_a) * (set_int + only_b) + (only_b + background) * (only_a + background)
  E <- E / background^2

  return(1 - ((O-E) / (1-E)))
}

#' Calculate the Kappa distance of all geneset combinations
#'
#' @param genesets a list of the genesets to score
#'
#' @return a [Matrix::Matrix()] with the pairwise Kappa distance of each geneset pair
#' @export
#'
#' @examples
#' genesets <- list(c("PDHB", "VARS2", "IARS2"), c("IARS2", "PDHA2"))
#' m <- getKappaMatrix(genesets)
getKappaMatrix <- function(genesets){
  l <- length(genesets)
  if(l == 0){
    return(-1)
  }
  k <- Matrix::Matrix(0, l, l)
  unique_genes <- length(unique(unlist(genesets)))

  for(i in 1:(l-1)){
    a <- genesets[i]
    for(j in (i+1):l){
      b <- genesets[j]
      k[i, j] <- k[j, i] <- calculateKappa(a, b, unique_genes)
    }
  }
  rownames(k) <- genesets$Geneset
  colnames(k) <- genesets$Geneset
  return(k)
}
