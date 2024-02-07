#' Get NCBI ID
#'
#' Get the NCBI ID of a species
#'
#' @param species character, the species of your input data
#' @param version character, the version of STRING you want to use, defaults to
#'                the current version of STRING
#'
#' @return A character of the NCBI ID of `species`
#' @export
#'
#' @examples
#' species <- "Homo sapiens"
#' id <- getId(species = species)
#'
#' species <- "Mus musculus"
#' id <- getId(species = species)
getId <- function(species, version = "11.5") {
  # Download available species information from STRING
  url_species <- sprintf("https://stringdb-static.org/download/species.v%s.txt",
                         version)

  # Read species data from URL
  df_species <- read.delim(url(url_species))

  # Get the species ID of the respective organism
  species_id <- df_species$X.taxon_id[match(species, df_species$official_name_NCBI)]

  # Return the species ID
  return(species_id)
}

#' Get the STRING db entry of a species
#'
#' Get the respective [STRINGdb] object of your species of interest
#'
#' @param species numeric, the NCBI ID of the species of interest
#' @param version character, The STRINGdb version to use, defaults to the
#'                current version
#' @param score_threshold numeric, A score threshold to cut the retrieved
#'                        interactions, defaults to 0 (all interactions)
#' @param cache_location character, the location to save the downloaded
#'                       [STRINGdb] object.
#'
#' @return a [STRINGdb] object of `species`
#' @export
#'
#' @import STRINGdb
#' @importFrom tools R_user_dir
#' @examples
#' species <- getId(species = "Homo sapiens")
#' string_db <- getStringDB(as.numeric(species))
getStringDB <- function(species,
                        version = "11.5",
                        score_threshold = 0.00,
                        cache_location = NULL) {

  if(is.null(cache_location)){
    cache_location <- tools::R_user_dir("GeDi", which = "cache")
    }

  return(STRINGdb$new(
        version = version,
        species = species,
        score_threshold = score_threshold,
        input_directory = cache_location
      ))
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
#' @param genes a `list`, A `list` of genes to download the respective protein-protein
#'              interaction information
#' @param string_db A [STRINGdb] object, the species of the object should match
#'                  the species of `genes`.
#' @param anno_df An annotation `data.frame` mapping [STRINGdb] ids to gene
#'                names, e.g. downloaded with \code{GeDi::getAnnotation()}
#'
#' @return A `data.frame` of Protein-Protein interactions
#' @export
#' @importFrom dplyr distinct
#' @import STRINGdb
#' @examples
#' ## Mock example showing how the data should look like
#'
#' genes <- c(c("CFTR", "RALA"), c("CACNG3", "ITGA3"), c("DVL2"))
#' string_db <- getStringDB(9606)
#' string_db
#' anno_df <- getAnnotation(string_db)
#' ppi <- getPPI(genes, string_db, anno_df)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#' string_db <- getStringDB(9606)
#' string_db
#' anno_df <- getAnnotation(string_db)
#' genes <- GeDi::getGenes(macrophage_topGO_example_small)
#' ppi <- getPPI(genes, string_db, anno_df)
getPPI <- function(genes, string_db, anno_df) {
  # Convert input list to vector
  genes <- unlist(genes)

  # Match gene names to STRINGdb IDs
  string_ids <- anno_df$STRING_id[match(genes, anno_df$alias)]

  # Get interaction scores of genes from the STRING database
  scores <- string_db$get_interactions(string_ids)
  colnames(scores) <- c("Gene1", "Gene2", "combined_score")

  # Normalize interaction scores to the (0, 1) interval
  max <- max(scores$combined_score,-Inf)
  min <- min(scores$combined_score, Inf)
  scores$combined_score <-
    round((scores$combined_score - min) / (max - min), 2)

  # Filter the data frame to ensure interactions are unique
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

  # Build up the final data frame of unique interactions
  scores <- distinct(scores)
  df <- distinct(df)

  # Return the final data frame of interactions and scores
  scores <- rbind(scores, df)
  return(scores)
}
