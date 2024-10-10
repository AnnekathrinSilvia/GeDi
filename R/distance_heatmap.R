#' Plot a heatmap
#'
#' Plot a heatmap of a matrix of (distance) scores of the input genesets
#'
#' @param distance_scores A [Matrix::Matrix()] of (distance) scores for each
#'                        pairwise combination of genesets.
#' @param chars_limit Numeric value, Indicates how many characters of the
#'                    row and column names of `distance_scores` should be
#'                    plotted. Defaults to 50 and prevents crowded axes due to
#'                    long names.
#' @param plot_labels Logical, Indicates if row and collabels should be plotted.
#'                    Defaults to TRUE
#' @param cluster_rows Logical, Indicates whether or not the rows should be 
#'                     clustered based on the distance scores. Defaults to TRUE
#' @param cluster_columns Logical, Indicates whether or not the rows should be 
#'                        clustered based on the distance scores. Defaults to 
#'                        TRUE
#' @param title character, a title for the figure. Defaults to "Distance Scores" 
#'              
#'
#' @return A [ComplexHeatmap::Heatmap()] plot object.
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @export
#'
#' @examples
#' ## Mock example showing how the data should look like
#'
#' distance_scores <- Matrix::Matrix(0.5, 20, 20)
#' distance_scores[c(11:15), c(2:6)] <- 0.2
#' rownames(distance_scores) <- colnames(distance_scores) <- as.character(c(1:20))
#' p <- distanceHeatmap(distance_scores)
#'
#' ## Example using the data available in the package
#' data(scores_macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' p <- distanceHeatmap(scores_macrophage_topGO_example_small)
distanceHeatmap <- function(distance_scores,
                            chars_limit = 50,
                            plot_labels = TRUE,
                            cluster_rows = TRUE,
                            cluster_columns = TRUE,
                            title = "Distance Scores") {
  # Check if distance scores are provided
  stopifnot(!is.null(distance_scores))
  stopifnot(chars_limit >= 0)

  if(plot_labels){
    # Cut the labels to the specified character limit
    labels <- substr(as.character(rownames(distance_scores)), 1, chars_limit)
    # Set truncated labels for row and column names
    rownames(distance_scores) <- colnames(distance_scores) <- labels
  }else{
    rownames(distance_scores) <- colnames(distance_scores) <- NULL
  }

  col_fun <- colorRamp2(c(0, 0.5, 1), c("red", "white", "blue"))
  
  # Create a heatmap using the distance scores matrix
  p <- Heatmap(as.matrix(distance_scores), 
               heatmap_legend_param = list(title = title),
               cluster_rows = cluster_rows,
               cluster_columns = cluster_columns)
  # Return the heatmap plot
  return(p)
}
