#' Split string of genes
#'
#' Split a long string of space separated genes into a `list` of individual
#' genes.
#'
#' @param genesets `data.frame`, a `data.frame` with at least two columns.
#'                 One should be called `Geneset`, containing the
#'                 names/identifiers of the genesets in the data. The second
#'                 column should be called `Genes` and contains one string of
#'                 the genes contained in each geneset.
#' @param gene_name character, alternative name for the column containing the
#'                  genes in `genesets`. If not given, the column is expected to
#'                  be called `Genes`.
#'
#' @return `list`, containing for each geneset in the `Geneset` column a
#'         `list` of the included genes.
#' @export
#'
#' @examples
#' df <- data.frame(
#'   Geneset = c(
#'     "Cell Cycle",
#'     "Biological Process",
#'     "Mitosis"
#'   ),
#'   Genes = c(
#'     c("PDHB,VARS2,IARS2"),
#'     c("LARS,LARS2"),
#'     c("IARS,SUV3")
#'   )
#' )
#' genes <- getGenes(df)
getGenes <- function(genesets, gene_name = NULL) {
  if (length(genesets) == 0) {
    return(NULL)
  }
  if (is.null(gene_name)) {
    stopifnot(any(names(genesets) == "Genes"))
  }

  # check in whcih column the genes are
  if (!is.null(gene_name)) {
    genesList <- genesets[, gene_name]
  } else {
    genesList <- genesets$Genes
  }

  # guess on the separator in the genes
  sep <- .findSeparator(genesList)
  # separate large string of genes into list of individual genes
  genes <- lapply(1:nrow(genesets), function(i) {
    toupper(strsplit(genesList[i], sep)[[1]])
  })

  return(genes)
}


#' Make an educated guess on the separator character
#'
#' This function tries to guess which separator was used in a list of delimited
#' strings.
#'
#' @param stringList `list`, a `list` of strings
#' @param sepList `list`, containing the candidates for being identified as
#'                 separators. Defaults to \code{c(",", "\t", ";"," ", "/")}.
#'
#' @return character, corresponding to the guessed separator. One of ","
#'         (comma), "\\t" (tab), ";" (semicolon)," " (whitespace) or "/"
#'         (backslash).
#'
#' @importFrom stringr str_count
.findSeparator <- function(stringList, sepList = c(",", "\t", ";", " ", "/")) {
  sephits_min <-
    sapply(sepList, function(x) {
      min(str_count(stringList, x))
    }) # minimal number of separators on all lines
  sep <- sepList[which.max(sephits_min)]

  return(sep)
}

#' Make an educated guess on the separator character
#'
#' This function tries to guess which separator was used in a text delimited
#' file.
#'
#' @param file character, location of a file to read data from.
#' @param sep_list `list`, containing the candidates for being identified as
#'                 separators. Defaults to \code{c(",", "\t", ";"," ", "/")}.
#'
#' @return character, corresponding to the guessed separator. One of ","
#'         (comma), "\\t" (tab), ";" (semicolon)," " (whitespace) or "/"
#'         (backslash).
.sepguesser <- function(file, sep_list = c(",", "\t", ";", " ", "/")) {
  rl <- readLines(file, warn = FALSE)
  rl <- rl[rl != ""] # allow last line to be empty
  sep <- .findSeparator(rl, sep_list)
  return(sep)
}

#' Filter Genesets from the input data
#'
#' Filter a preselected list of genesets from a `data.frame` of genesets
#'
#' @param remove `list`, a `list` of geneset identifiers to be removed from
#'               `df_genesets`
#' @param df_genesets `data.frame`, a `data.frame` of the genesets.
#'
#' @return A `data.frame` of the input data without the genesets listed in
#'         `remove`.
.filterGenesets <- function(remove,
                            df_genesets) {
  # get genesets to remove
  genesets_to_remove <- unlist(strsplit(remove, "\\s+"))
  results <- list()
  df_genesets <- df_genesets[!(df_genesets$Geneset %in% genesets_to_remove), ]

  # return information of the data.frame without the genesets in removed
  results[[1]] <- df_genesets
  results[[2]] <- df_genesets$Geneset
  genes <- getGenes(df_genesets)
  results[[3]] <- genes


  names(results) <- c("Geneset", "gs_names", "Genes")
  return(results)
}

#' Check if PPI has right format
#'
#' Check if the Protein-Protein-interaction (PPI) has the expected format for
#' this app
#'
#' @param ppi A PPI object
#'
#' @return A `data.frame` of the input `ppi` with the expected column names
#'
.checkPPI <- function(ppi) {
  stopifnot(is.data.frame(ppi))
  stopifnot(ncol(ppi) == 3)
  stopifnot(all(is.character(ppi[, 1])) & all(is.character(ppi[, 2])))
  stopifnot(all(is.numeric(ppi[, 3])))

  colnames(ppi) <- c("Gene1", "Gene2", "combined_score")
  return(ppi)
}

#' Check if iput genesets have the right format
#'
#' Check if the input genesets have the expected format of the app
#'
#' @param genesets A object of geneset input data
#'
.checkGenesets <- function(genesets) {
  stopifnot(is.data.frame(genesets))
  stopifnot(any(names(genesets) == "Genesets") & any(names(genesets) == "Genes"))
  stopifnot(all(is.character(genesets$Genesets)) & all(is.character(genesets$Genes)))
}

#' Determine the number of cores touse for a function
#'
#' Determine the number of cores to use for a function
#'
#' @param n_cores numeric, number of cores to use for the function.
#'                Defaults to `Null` in which case the function takes half of
#'                the available cores.
#'
#' @return The number of cores to use.
#'
.getNumberCores <- function(n_cores = NULL) {
  # check the number of cores to use
  available_cores <- parallel::detectCores()
  if (is.null(n_cores)) {
    n_cores <- max(round(available_cores/2), 1)
  } else {
    if (n_cores > available_cores) {
      n_cores <- available_cores - 1
    }
  }

  return(n_cores)
}

.actionButtonStyle <-
  "color: #FFFFFF; background-color: #0092AC; border-color: #0092AC"
.tourButtonStyle <-
  "color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"
