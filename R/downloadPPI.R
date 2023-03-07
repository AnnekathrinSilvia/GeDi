#' Get ID of a species
#'
#' Get the NCBI ID of a species
#'
#' @param species character, the species of your input data
#' @param version character, the version of STRING you want to use, defaults to
#'                11.5, the current version of STRING
#'
#' @return a character of the NCBI ID of `species`
#' @export
#'
#' @examples
#' species <- "Homo sapiens"
#' id <- getId(species = species)
getId <- function(species, version = "11.5") {
  # download available species from STRING
  url_species <- sprintf(
    "https://stringdb-static.org/download/species.v%s.txt",
    version
  )
  df_species <- read.delim(url(url_species))

  # get species ID of respective organism
  species_id <- df_species$X.taxon_id[
    match(species, df_species$official_name_NCBI)
  ]

  return(species_id)
}


#' Get the STRING entry of a species
#'
#' Get the respective [STRINGdb] object of your species of interest
#'
#' @param species numeric, the NCBI ID of the species of interest
#' @param version character, the STRINGdb version to use, defaults to 11.5
#'                (current version)
#' @param score_threshold numeric, a score threshold to cut the retrieved
#'                        interactions, defaults to 0 (all interactions)
#'
#' @return a [STRINGdb] object `species`
#' @export
#'
#' @import STRINGdb
#' @examples
#' species <- getId(species = "Homo Sapiens")
#' string_db <- getStringDB(as.numeric(species))
getStringDB <- function(species,
                        version = "11.5",
                        score_threshold = 0.00) {
  return(
    STRINGdb$new(
      version = version,
      species = species,
      score_threshold = score_threshold,
      input_directory = ""
    )
  )
}

#' Get the annotation of a [STRINGdb] object
#'
#' Get the annotation of a [STRINGdb] object, i.e. the aliases of the protein
#' information
#'
#' @param stringdb the [STRINGdb] object
#'
#' @return A `data.frame` mapping [STRINGdb] ids to gene names
#' @export
#'
#' @import STRINGdb
#' @examples
#' string_db <- getStringDB(9606)
#' string_db
#' anno_df <- getAnnotation(string_db)
getAnnotation <- function(stringdb) {
  return(stringdb$get_aliases())
}

#' Download Protein-Protein Interaction (PPI)
#'
#' Download the Protein-Protein Interaction (PPI) information of a [STRINGdb]
#' object
#'
#' @param genes list, a list of genes to download the respective protein-protein
#'              interaction information
#' @param string_db A [STRINGdb] object, the species of the object should match
#'                  the species of `genes`
#' @param anno_df A annotation `data.frame` mapping [STRINGdb] ids to gene
#'                names, e.g. downloaded with [GeDi::getAnnotation()]
#'
#' @return A `data.frame` of Protein-Protein interactions
#' @export
#' @importFrom dplyr distinct
#' @import STRINGdb
#' @examples
#'
#' genes <- c(c("CFTR", "RALA"), c("CACNG3", "ITGA3"), c("DVL2"))
#' string_db <- getStringDB(9606)
#' string_db
#' anno_df <- getAnnotation(string_db)
#' ppi <- getPPI(genes, string_db, anno_df)
getPPI <- function(genes, string_db, anno_df) {
  genes <- unlist(genes)

  # match gene names to STRINGdb ids
  string_ids <- anno_df$STRING_id[match(genes, anno_df$alias)]
  # get interaction scores of genes
  scores <- string_db$get_interactions(string_ids)
  colnames(scores) <- c("Gene1", "Gene2", "combined_score")

  max <- max(scores$combined_score, -Inf)
  min <- min(scores$combined_score, Inf)

  # normalize the interaction scores to (0, 1)
  scores$combined_score <-
    round((scores$combined_score - min) / (max - min), 2)

  # filter the final data frame as interactions are symmetric
  gene_names_to <-
    anno_df$alias[match(scores$Gene2, anno_df$STRING_id)]
  gene_names_from <-
    anno_df$alias[match(scores$Gene1, anno_df$STRING_id)]

  scores$Gene2 <- gene_names_to
  scores$Gene1 <- gene_names_from

  reverse_to <- scores$Gene1
  reverse_from <- scores$Gene2

  df <-
    data.frame(
      Gene1 = reverse_from,
      Gene2 = reverse_to,
      combined_score = scores$combined_score
    )

  # build up final data frame of unique interactions
  scores <- distinct(scores)
  df <- distinct(df)

  scores <- rbind(scores, df)
  return(scores)
}
