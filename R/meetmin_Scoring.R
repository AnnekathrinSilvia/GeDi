#' Calculate the Meet-Min distance of all geneset combinations
#'
#' @param genesets a list of the genesets to score
#'
#' @return a [Matrix::Matrix()] with the pairwise Meet-Min distance of each geneset pair
#' @export
#'
#' @examples
#' genesets <- list(c("PDHB", "VARS2", "IARS2"), c("IARS2", "PDHA2"))
#' m <- getMeetMinMatrix(genesets)
getMeetMinMatrix <- function(genesets){
  l <- length(genesets)
  if(l == 0){
    return(-1)
  }
  m <- Matrix::Matrix(0, l, l)

  for(i in 1:(l-1)){
    a <- genesets[i]
    for(j in (i+1):l){
      b <- genesets[j]
      int <- length(intersect(a, b))

      m[i, j] <- m[j, i] <- 1 - (int / min(length(a), length(b)))
    }
  }

  rownames(m) <- genesets$Geneset
  colnames(m) <- genesets$Geneset
  return(m)
}
