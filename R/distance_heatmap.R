#' Plot a heatmap on the distance scores of the genesets
#'
#' Plot a heatmap of a matrix of (distance) scores of the input genesets
#'
#' @param distance_scores A [Matrix::Matrix()] of (distance) scores for each
#'                        pairwise combination of genesets.
#' @param chars_limit Numeric value, how many characters of the geneset names
#'                    of `distance_scores`. Defaults to 50.
#'
#' @return A plot returned by the [ggplot2::ggplot()] function
#' @import ggplot2
#' @import viridis
#' @export
#'
#' @examples
#' distance_scores <- Matrix::Matrix(0.5, 20, 20)
#' distance_scores[c(11:15), c(2:6)] <- 0.2
#' rownames(distance_scores) <- colnames(distance_scores) <- as.character(c(1:20))
#' p <- distance_heatmap(distance_score, cluster = T)
distance_heatmap <- function(distance_scores,
                             chars_limit = 50,
                             hcluster = FALSE){
  stopifnot(!is.null(distance_scores))

  labels <- substr(as.character(rownames(distance_scores)), 1, chars_limit)

  df <- expand.grid(Geneset1 = rownames(distance_scores), Geneset2 = colnames(distance_scores))
  scores <- unlist(as.list(as.matrix(distance_scores)))
  df$distance_score <- scores

  if(hcluster){
    m <- tidyr::pivot_wider(df, names_from = "Geneset1", values_from = "distance_score")
    m <- as.matrix(m[, -1])
    ord <- hclust(dist(t(m)), method = "ward.D")$order
  }

  p <- ggplot(df, aes(Geneset1, Geneset2, fill = distance_score)) +
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    hrbrthemes::theme_ipsum() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

  if(hcluster){
   p <- p + scale_y_discrete(limits = colnames(m)[ord])
  }

  return(p)
}
