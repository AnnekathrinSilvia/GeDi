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
distanceHeatmap_OLD <- function(distance_scores,
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
#' @param display_similarity Logical, Indicates whether or not the scores should be 
#'                           plotted as a distance matrix (i.e. 0 indicates completely
#'                           identical sets) or as a similarity matrix (i.e. 1
#'                           indicates completely identical sets). Defaults to FALSE
#'                           (i.e. a distance matrix).
#'                    
#' @param quantile_limits Numerical vector, Used to scale the colors in the 
#'                        heatmap between the given quantiles. Can be helpful in
#'                        order to reduce the effect of the diagonal or other 
#'                        outlier values on the color scheme. Defaults to 
#'                        c(0.025, 0.975).
#'              
#'
#' @return A [ComplexHeatmap::Heatmap()] plot object.
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom stats quantile
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
                            title = "Distance Scores",
                            display_similarity = FALSE,
                            quantile_limits = c(0.025, 0.975)) {
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
  
  # Set the color function for the matrix, also scale it only to the given 
  # qauntiles
  q_min <- quantile(distance_scores, quantile_limits[1], na.rm = TRUE)
  q_max <- quantile(distance_scores, quantile_limits[2], na.rm = TRUE)
  col_mid <- colorRamp2(
    c(q_min, (q_min+q_max)/2, q_max),
    c("white", "lightblue", "darkblue")
  )
  # Set the final color function
  col_fun <- function(x) {
    x[x < q_min] <- q_min
    x[x > q_max] <- q_max 
    col_mid(x)
  }
  
  if(display_similarity){
    # Create a heatmap using the similarity scores matrix
    similarity_scores <- as.matrix(1-distance_scores)
    
    # Set the color function for the matrix, also scale it only to the given 
    # qauntiles
    q_min <- quantile(similarity_scores, quantile_limits[1], na.rm = TRUE)
    q_max <- quantile(similarity_scores, quantile_limits[2], na.rm = TRUE)
    col_mid <- colorRamp2(
      c(q_min, (q_min+q_max)/2, q_max),
      c("white", "red", "darkred")
    )
    # Set the final color function
    col_fun <- function(x) {
      x[x < q_min] <- q_min 
      x[x > q_max] <- q_max  
      col_mid(x)
    }
    
    p <- Heatmap(as.matrix(similarity_scores), 
                 heatmap_legend_param = list(title = "Similarity Scores"),
                 cluster_rows = cluster_rows,
                 cluster_columns = cluster_columns,
                 col = col_fun)
  }else{
    # Create a heatmap using the distance scores matrix
    p <- Heatmap(as.matrix(distance_scores), 
                 heatmap_legend_param = list(title = title),
                 cluster_rows = cluster_rows,
                 cluster_columns = cluster_columns,
                 col = col_fun)
  }
  # Return the heatmap plot
  return(p)
}
