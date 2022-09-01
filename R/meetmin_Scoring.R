#' Get Matrix of Meet-Min distances
#'
#' Calculate the Meet-Min distance of all combinations of genesets in a given
#' data set of genesets.
#'
#' @param genesets A `list` of genesets (each geneset is represented by a `list`
#'                 of the corresponding genes).
#' @param progress An optional [shiny::Progress()] object to track the progress
#'                 progress of the function in the app.
#'
#' @return A [Matrix::Matrix()] with the pairwise Meet-Min distance of each
#'         geneset pair.
#' @export
#'
#' @examples
#' genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
#' m <- getMeetMinMatrix(genesets)
getMeetMinMatrix <- function(genesets, progress = NULL){
  l <- length(genesets)
  if(l == 0){
    return(-1)
  }
  m <- Matrix::Matrix(0, l, l)

  for(i in 1:(l - 1)){
    a <- genesets[i]
    for(j in (i+1):l){
      b <- genesets[j]

      if(length(a) == 0 || length(b) == 0){
        m[i, j] <- m[j, i] <- 1
      }else{
        int <- length(intersect(a, b))
        m[i, j] <- m[j, i] <- 1 - (int / min(length(a), length(b)))
      }
    }
    if(!is.null(progress)){
      progress$inc(1/l, detail = paste("Scoring geneset number", i))
    }
  }
  return(m)
}
