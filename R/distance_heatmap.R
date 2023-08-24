#' Plot a heatmap
#'
#' Plot a heatmap of a matrix of (distance) scores of the input genesets
#'
#' @param distance_scores A [Matrix::Matrix()] of (distance) scores for each
#'                        pairwise combination of genesets.
#' @param chars_limit Numeric value, Indicates how many characters of the
#'                    row and column names of `distance_scores` should be
#'                    plotted. Defaults to 50 and prevents crowded axises due to
#'                    long names.
#'
#' @return A [ComplexHeatmap::Heatmap()] plot object.
#' @importFrom ComplexHeatmap Heatmap
#' @export
#'
#' @examples
#' distance_scores <- Matrix::Matrix(0.5, 20, 20)
#' distance_scores[c(11:15), c(2:6)] <- 0.2
#' rownames(distance_scores) <- colnames(distance_scores) <- as.character(c(1:20))
#' p <- distance_heatmap(distance_scores)
distance_heatmap <- function(distance_scores, chars_limit = 50) {
  # Check if distance scores are provided
  stopifnot(!is.null(distance_scores))

  # Cut the labels to the specified character limit
  labels <- substr(as.character(rownames(distance_scores)), 1, chars_limit)

  # Set truncated labels for row and column names
  rownames(distance_scores) <- colnames(distance_scores) <- labels

  # Create a heatmap using the distance scores matrix
  p <- Heatmap(as.matrix(distance_scores),
               name = "Distance Scores"
  )

  # Return the heatmap plot
  return(p)
}
