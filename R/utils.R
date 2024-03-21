#' Split string of genes
#'
#' Split a long string of space separated genes into a `list` of individual
#' genes.
#'
#' @param genesets a `data.frame`, A `data.frame` with at least two columns.
#'                 One should be called `Geneset`, containing the
#'                 names/identifiers of the genesets in the data. The second
#'                 column should be called `Genes` and contains one string of
#'                 the genes contained in each geneset.
#' @param gene_name a character, Alternative name for the column containing the
#'                  genes in `genesets`. If not given, the column is expected to
#'                  be called `Genes`.
#'
#' @return A `list` containing for each geneset in the `Geneset` column a
#'         `list` of the included genes.
#' @export
#'
#' @examples
#' ## Mock example showing how the data should look like
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
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' genes <- getGenes(macrophage_topGO_example_small)
getGenes <- function(genesets, gene_name = NULL) {
  # If there are no genesets, return NULL
  if (length(genesets) == 0) {
    return(NULL)
  }

  # If gene_name is not provided, ensure that a "Genes" column exists
  if (is.null(gene_name)) {
    stopifnot(any(names(genesets) == "Genes"))
  }

  # Determine in which column the gene information is stored
  if (!is.null(gene_name)) {
    genesList <- genesets[, gene_name]
  } else {
    genesList <- genesets$Genes
  }

  # Guess the separator used in the gene lists
  sep <- .findSeparator(genesList)

  # Split large strings of genes into individual gene lists
  genes <- lapply(seq_len(nrow(genesets)), function(i) {
    toupper(strsplit(genesList[i], sep)[[1]])
  })

  # Return the list of extracted gene sets
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
#' @references
#' See https://github.com/federicomarini/ideal for details on the original
#' implementation.
#'
#' @return character, corresponding to the guessed separator. One of ","
#'         (comma), "\\t" (tab), ";" (semicolon)," " (whitespace) or "/"
#'         (backslash).
#'
.findSeparator <- function(stringList, sepList = c(",", "\t", ";", " ", "/")) {
  sephits_min <-
    vapply(sepList, function(x) {
      sum(vapply(stringList, function(y) nchar(y) - nchar(gsub(x, '', y)), numeric(1)))
    }, numeric(1))
  sep <- sepList[which.max(sephits_min)]

  return(sep)
}

#' Make an educated guess on the separator character
#'
#' This function tries to guess which separator was used in a text delimited
#' file.
#'
#' @param file a character, location of a file to read data from.
#' @param sep_list a `list`, containing the candidates for being identified as
#'                 separators. Defaults to \code{c(",", "\t", ";"," ", "/")}.
#'
#' @references
#' See https://github.com/federicomarini/ideal for details on the original
#' implementation.
#'
#' @return A character, corresponding to the guessed separator. One of ","
#'         (comma), "\\t" (tab), ";" (semicolon)," " (whitespace) or "/"
#'         (backslash).
.sepguesser <- function(file, sep_list = c(",", "\t", ";", " ", "/")) {
  rl <- readLines(file, warn = FALSE)
  # Allow last line to be empty
  rl <- rl[rl != ""]
  sep <- .findSeparator(rl, sep_list)
  return(sep)
}

#' Filter Genesets from the input data
#'
#' Filter a preselected list of genesets from a `data.frame` of genesets
#'
#' @param remove a `list`, A list of geneset names to be removed
#' @param df_genesets a `data.frame`, A `data.frame` with at least two columns.
#'                 One should be called `Geneset`, containing the
#'                 names/identifiers of the genesets in the data. The second
#'                 column should be called `Genes` and contains one string of
#'                 the genes contained in each geneset.
#'
#' @return A `data.frame` containing information about filtered genesets
.filterGenesets <- function(remove,
                            df_genesets) {
  # Get genesets to remove
  if (length(remove) > 0) {
    # Split the remove vector
    genesets_to_remove <- unlist(remove)
    df_genesets <- df_genesets[!(df_genesets$Geneset %in% genesets_to_remove), ]
  }

  # Initialize a list to store results
  results <- list()
  # Store the data frame without the removed genesets
  results[[1]] <- df_genesets

  # Store the names of the remaining genesets
  results[[2]] <- df_genesets$Geneset

  # Extract gene information for the remaining genesets
  genes <- getGenes(df_genesets)
  results[[3]] <- genes

  # Rename the elements in the results list
  names(results) <- c("Geneset", "gs_names", "Genes")

  # Return the filtered geneset information
  return(results)
}

#' Title
#'
#' @param genesets a `data.frame`, A `data.frame` with at least two columns.
#'                 One should be called `Geneset`, containing the
#'                 names/identifiers of the genesets in the data. The second
#'                 column should be called `Genes` and contains one string of
#'                 the genes contained in each geneset.
#'
#' @return a `list` of geneset descriptions
.getGenesetDescriptions <- function(genesets){
  columns <- names(genesets)
  terms <- list()
  if ("Term" %in% columns) {
    terms <- genesets$Term
  } else if ("Description" %in% columns) {
    terms <- genesets$Description
  } else {
    terms <- genesets$Genesets
  }

  if(any(duplicated(terms))){
    duplicates_ids <- which(duplicated(terms))
    duplicates <- terms[duplicates_ids]
    duplicates <- paste(duplicates, "1", sep = "_")
    terms[duplicates_ids] <- duplicates
  }

  return(terms)
}

#' Check PPI format
#'
#' Check if the Protein-Protein-interaction (PPI) has the expected format for
#' this app
#'
#' @param ppi a `data.frame`, Protein-protein interaction (PPI) network data frame.
#'            The object is expected to have three columns, `Gene1` and `Gene2`
#'            which specify the gene names of the interacting proteins in no
#'            particular order (symmetric interaction) and a column
#'            `combined_score` which is a numerical value of the strength of
#'            the interaction.
#'
#' @return A validated and formatted PPI data frame.
#'
.checkPPI <- function(ppi) {
  # Check if ppi is a data frame
  stopifnot(is.data.frame(ppi))

  # Check if ppi has the correct number of columns
  stopifnot(ncol(ppi) == 3)

  # Check if Gene1 and Gene2 columns are character type
  stopifnot(all(is.character(ppi[, 1])) & all(is.character(ppi[, 2])))

  # Check if combined_score column is numeric
  stopifnot(all(is.numeric(ppi[, 3])))

  # Rename columns to expected names
  colnames(ppi) <- c("Gene1", "Gene2", "combined_score")

  # Return the validated and formatted PPI data frame
  return(ppi)
}

#' Check genesets format
#'
#' Check if the input genesets have the expected format for this app
#'
#' @param genesets a `list`, A `list` of genesets where each genesets is represented
#'                 by `list` of genes.
#'
#' @return A validated and formatted genesets data frame.
#'
.checkGenesets <- function(genesets) {
  # Check if genesets is a data frame
  stopifnot(is.data.frame(genesets))

  # Check if genesets has the required columns
  stopifnot(any(names(genesets) == "Genesets") & any(names(genesets) == "Genes"))

  # Check if Genesets and Genes columns are character type
  stopifnot(all(is.character(genesets$Genesets)) & all(is.character(unlist(genesets$Genes))))

  # Return the validated and formatted genesets data frame
  return(genesets)
}

#' Check distance scores format
#'
#' Check if the provided distance scores have the expected format for this app
#'
#' @param genesets a `list`, A `list` of genesets where each genesets is represented
#'                 by `list` of genes.
#' @param distance_scores A [Matrix::Matrix()] or  object,
#'                        A matrix with numerical (distance) scores.
#'
#' @return A validated and formatted distance_scores [Matrix::Matrix()].
#'
.checkScores <- function(genesets, distance_scores){
  # Check if the distance_scores matrix is square
  stopifnot(nrow(distance_scores) == ncol(distance_scores))

  # Check if the number of rows in genesets matches the number of rows
  # in distance_scores
  stopifnot(nrow(genesets) == nrow(distance_scores))

  # Check if all values in the distance_scores matrix are numeric
  stopifnot(all(is.numeric(distance_scores@x)))

  # Check if row names and column names in distance_scores are the same
  stopifnot(rownames(distance_scores) == colnames(distance_scores))
  # Check if row names in distance_scores match the Genesets column in the
  # genesets data frame
  stopifnot(rownames(distance_scores) == genesets$Genesets)

  # Return the validated and formatted distance_scores matrix
  return(distance_scores)
}

#' Determine the number of cores to use for a function
#'
#' Determine the number of CPU cores the scoring functions should use when
#' computing the distance scores.
#'
#' @param n_cores numeric, number of cores to use for the function.
#'                Defaults to `Null` in which case the function takes half of
#'                the available cores.
#'
#' @return Number of CPU cores to be used.
.getNumberCores <- function(n_cores = NULL) {
  # Check the number of available CPU cores
  available_cores <- parallel::detectCores()

  # If n_cores is not provided, set it to half of available_cores (minimum 1)
  if (is.null(n_cores)) {
    n_cores <- max(round(available_cores / 2), 1)
  } else {
    # If n_cores exceeds available_cores, adjust n_cores to available_cores - 1
    if (n_cores > available_cores) {
      n_cores <- available_cores - 1
    }
  }

  # Return the determined number of CPU cores
  return(n_cores)
}

# Shiny resource paths ----------------------------------------------------

.onLoad <- function(libname, pkgname) {
  # Create link to logo
  shiny::addResourcePath("GeDi", system.file("www", package = "GeDi"))
}

.actionButtonStyle <-
  "color: #FFFFFF; background-color: #0092AC; border-color: #0092AC"
#.tourButtonStyle <- "color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"

