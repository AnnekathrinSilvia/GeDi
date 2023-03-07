#' Plot a dendrogram on the distance scores of the genesets
#'
#' Plot a dendrogram of a matrix of (distance) scores of the input genesets
#'
#' @param distance_scores A [Matrix::Matrix()] containing (distance) scores
#'                        between 0 and 1.
#' @param cluster_method character, indicating the clustering method
#'                       for the [stats::hclust()] function. See the
#'                       [stats::hclust()] function for the available options.
#'                       Defaults to 'average'.
#'
#' @return A dendrogram returned by the [ggdendro::ggdendrogram()] function.
#' @export
#' @importFrom ggdendro ggdendrogram
#' @importFrom stats dist hclust
#'
#' @examples
#' distance_scores <- Matrix::Matrix(0.5, 20, 20)
#' distance_scores[c(11:15), c(2:6)] <- 0.2
#' dendro <- distance_dendro(distance_scores, cluster_method = "single")
distance_dendro <- function(distance_scores,
                            cluster_method = "average") {
  dist <- dist(t(distance_scores))
  hc <- hclust(dist, cluster_method)
  p <- ggdendrogram(hc, rotate = FALSE, size = 2)

  return(p)
}
