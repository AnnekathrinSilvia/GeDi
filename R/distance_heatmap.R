#' Plot a heatmap on the distance scores of the genesets
#'
#' Plot a heatmap of a matrix of (distance) scores of the input genesets
#'
#' @param distance_scores A [Matrix::Matrix()] of (distance) scores for each
#'                        pairwise combination of genesets.
#' @param chars_limit Numeric value, how many characters of the geneset names
#'                    of `distance_scores`. Defaults to 50.
#'
#' @return A plot returned by the [ComplexHeatmap::Heatmap()] function
#' @importFrom ComplexHeatmap Heatmap
#' @export
#'
#' @examples
#' distance_scores <- Matrix::Matrix(0.5, 20, 20)
#' distance_scores[c(11:15), c(2:6)] <- 0.2
#' rownames(distance_scores) <- colnames(distance_scores) <- as.character(c(1:20))
#' p <- distance_heatmap(distance_scores)
distance_heatmap <- function(distance_scores,
                             chars_limit = 50) {
  # check if there are distance scores
  stopifnot(!is.null(distance_scores))

  # cut the labels to chars_limit
  labels <- substr(as.character(rownames(distance_scores)), 1, chars_limit)

  # set cut labels for plot annotation
  rownames(distance_scores) <- colnames(distance_scores) <- labels
  p <- Heatmap(as.matrix(distance_scores),
    name = "Distance Scores"
  )

  return(p)
}

