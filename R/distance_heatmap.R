#' Plot a heatmap on the distance scores of the genesets
#'
#' Plot a heatmap of a matrix of (distance) scores of the input genesets
#'
#' @param distance_scores A [Matrix::matrix()] of (distance) scores for each
#'                        pairwise combination of genesets.
#' @param chars_limit Numeric value, how many characters of the geneset names
#'                    of `distance_scores`. Defaults to 50.
#'
#' @return A plot returned by the [pheatmap::pheatmap()] function
#' @export
#'
#' @examples
distance_heatmap <- function(distance_scores, chars_limit = 50){
  stopifnot(!is.null(distance_scores))

  labels <- substr(as.character(rownames(distance_scores)), 1, chars_limit)

  p <- pheatmap::pheatmap(distance_scores,
                          color = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100),
                          labels_row = labels,
                          labels_col = labels)

  return(p)
}
