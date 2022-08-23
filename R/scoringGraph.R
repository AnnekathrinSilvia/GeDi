#' Title
#'
#' @param distanceMatrix a Distance Matrix containing for each pair the distance
#' @param cutOff a cutOff value classifying a threshold for similarity
#'
#' @return
#' @import Matrix
#' @export
#'
#' @examples
getAdjacencyMatrix <- function(distanceMatrix, cutOff){
  l <- nrow(distanceMatrix)
  adjMat <- Matrix::Matrix(0, l, l)

  for(i in 1:l){
    edge <- which(distanceMatrix[i, ] <= cutOff)
    no_edge <- which(distanceMatrix[i, ] > cutOff)

    adjMat[i, edge] <- 1
    adjMat[i, no_edge] <- 0
  }
  rownames(adjMat) <- rownames(distanceMatrix)
  colnames(adjMat) <- colnames(distanceMatrix)
  return(adjMat)
}

#' Title
#'
#' @param adjMatrix
#'
#' @return
#' @importFrom igraph graph_from_adjacency_matrix
#' @export
#'
#' @examples
buildGraph <- function(adjMatrix){
  g <- igraph::graph_from_adjacency_matrix(
    adjMatrix,
    mode = "undirected",
    add.colnames = NULL,
    add.rownames = NA
  )
}
