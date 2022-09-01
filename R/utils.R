#' Split string of genes
#'
#' Split a long string of space separated genes into a `list` of individual
#' genes.
#'
#' @param genesets A `data.frame` with at least two columns. One should be
#'                 called `Geneset`, which contains the names of the genesets
#'                 in the data (or any other interpretative identifier for the
#'                 genesets). The second column should be called `Genes` and
#'                 contains one string of the genes contained in each geneset.
#'
#' @return A `list` which contains for each geneset in the `Geneset` column a
#'         `list` of the included genes.
#' @export
#'
#' @examples
#' df <- data.frame(Geneset = c("Cell Cycle",
#'                            "Biological Process",
#'                             "Mitosis"),
#'                  Genes =  c(c("PDHB,VARS2,IARS2"),
#'                           c("LARS,LARS2"),
#'                           c("IARS,SUV3")))
#' genes <- getGenes(df)
getGenes <- function(genesets) {
  if (length(genesets) == 0) {
    return(list())
  }
  stopifnot(any(names(genesets) == "Genes"))

  genes <- lapply(1:nrow(genesets), function(i) {
    toupper(strsplit(genesets$Genes[i], ",")[[1]])
  })

  return(genes)
}


#' Make an educated guess on the separator character
#'
#' This function tries to guess which separator was used in a text delimited
#' file.
#'
#' @param file The name of the file which the data are to be read from.
#' @param sep_list A vector containing the candidates for being identified as
#'                 separators. Defaults to \code{c(",", "\t", ";"," ")}.
#'
#' @return A character value, corresponding to the guessed separator. One of ","
#'         (comma), "\\t" (tab), ";" (semicolon)," " (whitespace).
#' @export
#' @importFrom stringr str_count
#'
#' @examples
#' sepguesser(system.file("extdata/design_commas.txt", package = "ideal"))
#' sepguesser(system.file("extdata/design_semicolons.txt", package = "ideal"))
#' sepguesser(system.file("extdata/design_spaces.txt", package = "ideal"))
#' mysep <- sepguesser(system.file("extdata/design_tabs.txt",
#'                     package = "ideal"))
#'
sepguesser <- function(file, sep_list = c(",", "\t", ";", " ")) {
  separators_list <- sep_list
  rl <- readLines(file, warn = FALSE)
  rl <- rl[rl != ""] # allow last line to be empty
  sephits_min <-
    sapply(separators_list, function(x)
      min(stringr::str_count(rl, x))) # minimal number of separators on all lines
  sep <- separators_list[which.max(sephits_min)]
  return(sep)
}


.actionButtonStyle <- "color: #FFFFFF; background-color: #0092AC; border-color: #0092AC"
.tourButtonStyle <- "color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"
