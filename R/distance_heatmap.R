#' Title
#'
#' @param distance_scores
#' @param chars_limit
#'
#' @return
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
