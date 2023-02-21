#' Title
#'
#' @param distance_scores
#' @param cluster_method
#'
#' @return
#' @export
#' @importFrom ggdendro ggdendrogram
#'
#' @examples
distance_dendro <- function(distance_scores,
                            cluster_method = "average"){
  dist <- dist(t(distance_scores))
  hc <- hclust(dist, cluster_method)
  p <- ggdendrogram(hc, rotate = FALSE, size = 2)

  return(p)
}
