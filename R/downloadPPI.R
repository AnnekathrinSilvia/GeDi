#' Get NCBI ID of a species
#'
#' @param species the species you want to get the NCBI ID of
#'
#' @return the NCBI ID of your searched species
#' @export
#' @importFrom taxize get_ids
#'
#' @examples
#' id <- getId(species = "Homo Sapiens")
# getId <- function(species) {
#   id <- get_ids(species, db = "ncbi")
#   return(id$ncbi[species])
# }


#' Title
#'
#' @param version
#'
#' @return
#' @export
#'
#' @examples
getId <- function(species, version = "11.5") {
  # see from https://string-db.org/cgi/download

  ## accessory data
  ## e.g. https://stringdb-static.org/download/species.v11.5.txt

  url_species <- sprintf(
    "https://stringdb-static.org/download/species.v%s.txt",
    version
  )

  df_species <- read.delim(url(url_species))

  species_id <- df_species$X.taxon_id[
    match(species, df_species$official_name_NCBI)]

  return(species_id)


}


#' Get the STRINGdb entry of your searched species
#'
#' @param species the NCBI ID of the species you want to search (numeric)
#' @param version the STRINGdb version you want to use (defaults to the current version) (string)
#' @param score_threshold a score threshold to cut out interactions with a lower score (defaults to 0/all interactions) (numeric)
#' @param input_directory the directory to save the `STRINGdb` object
#'
#' @return a `STRINGdb` object of your searched species
#' @export
#'
#' @import STRINGdb
#' @examples
#' species <- getId(species = "Homo Sapiens")
#' string_db <- getStringDB(as.numeric(species))
getStringDB <- function(species, version = "11.5", score_threshold = 0.00, input_directory = "") {
  return(
    STRINGdb$new(
      version = version,
      species = species,
      score_threshold = score_threshold,
      input_directory = input_directory
    )
  )
}

#' Get the annotation of a `STRINGdb` object
#'
#' @param stringdb the `STRINGdb` object to get the annotation from
#'
#' @return a `data.frame` which maps `STRING_ids` to `alias`, which are gene names
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

#' Download the Protein-Protein Interaction information form you `STRINGdb` object
#'
#' @param genes the genes which are interested in
#' @param string_db a `STRINGdb` onject of the species you genes belong to
#' @param anno_df a annotation `data.frame` which maps `STRING_ids` to gene names
#'
#' @return a `data.frame` of Protein-Protein interactions
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
  string_ids <- anno_df$STRING_id[match(genes, anno_df$alias)]
  scores <- string_db$get_interactions(string_ids)
  colnames(scores) <- c("Gene1", "Gene2", "combined_score")

  max <- max(scores$combined_score, -Inf)
  min <- min(scores$combined_score, Inf)

  scores$combined_score <-
    round((scores$combined_score - min) / (max - min), 2)

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
  scores <- distinct(scores)
  df <- distinct(df)

  scores <- rbind(scores, df)
  return(scores)
}
