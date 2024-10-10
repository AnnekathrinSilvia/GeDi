#' Deprecated functions in GeDi
#'
#' Functions that are on their way to the function afterlife.
#' Their successors are also listed.
#'
#' The successors of these functions are likely coming from a renaming of the
#' functions to more intuitive function names
#'
#' @param ... Ignored arguments.
#'
#' @return All functions throw a warning, with a deprecation message pointing
#' towards its descendent (if available).
#'
#' @name deprecated
#'
#' @section Renaming funciton with more intuitive names:
#'
#' - [getGenes()], now replaced by the more intuitive name
#' [prepareGenesetData()]. The only change in its functionality concerns the
#' function name.
#'
#' @author Annekathrin Nedwed
#'
#' @examples
#' # try(getGenes())
#'
NULL
## #' @export
## #' @rdname defunct
## trendVar <- function(...) {
##   .Defunct("fitTrendVar")
## }
