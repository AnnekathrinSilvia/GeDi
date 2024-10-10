#' Plot a dendrogram
#'
#' Plot a dendrogram of a matrix of (distance) scores.
#'
#' @param distance_scores A [Matrix::Matrix()] containing (distance) scores
#'                        between 0 and 1.
#' @param cluster_method character, indicating the clustering method
#'                       for the [stats::hclust()] function. See the
#'                       [stats::hclust()] function for the available options.
#'                       Defaults to 'average'.
#'
#' @return A [ggdendro::ggdendrogram()] plot object.
#' @export
#' @importFrom ggdendro ggdendrogram
#' @importFrom stats dist hclust
#' @import Matrix
#'
#' @examples
#'
#' ## Mock example showing how the data should look like
#'
#' distance_scores <- Matrix::Matrix(0.5, 20, 20)
#' distance_scores[c(11:15), c(2:6)] <- 0.2
#' dendro <- distanceDendro(distance_scores, cluster_method = "single")
#'
#' ## Example using the data available in the package
#' data(scores_macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' dendro <- distanceDendro(scores_macrophage_topGO_example_small,
#'                         cluster_method = "average")
distanceDendro <- function(distance_scores,
                           cluster_method = "average") {
  # Check if distance scores are provided
  stopifnot(length(distance_scores) > 0 &&
              (!is.null(distance_scores)))
  
  # Calculate the distance between the distance scores and perform
  # hierarchical clustering
  dist <- dist(t(distance_scores))
  hc <- hclust(dist, cluster_method)
  # Create and plot the dendrogram
  p <- ggdendrogram(hc, rotate = FALSE, size = 2)
  # Return the dendrogram plot
  return(p)
}
