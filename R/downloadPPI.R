#' Title
#'
#' @param species
#'
#' @return
#' @export
#'
#' @examples
getId <- function(species) {
  id <- taxize::get_ids(species, db = "ncbi")
  return(id$ncbi[species])
}

#' Title
#'
#' @param species
#' @param version
#' @param score_threshold
#'
#' @return
#' @export
#'
#' @examples
getStringDB <- function(species, version, score_threshold) {
  return(
    STRINGdb$new(
      version = version,
      species = species,
      score_threshold = score_threshold,
      input_directory = ""
    )
  )
}

#' Title
#'
#' @param stringdb
#'
#' @return
#' @export
#'
#' @examples
getAnnotation <- function(stringdb) {
  return(stringdb$get_aliases())
}

#' Title
#'
#' @param genes
#' @param string_db
#' @param anno_df
#'
#' @return
#' @export
#'
#' @examples
getPPI <- function(genes, string_db, anno_df) {
  l <- length(genes)
  string_ids <- anno_df$STRING_id[match(genes, anno_df$alias)]
  scores <- string_db$get_interactions(string_ids)
  max <- max(scores$combined_score)
  min <- min(scores$combined_score)

  scores$combined_score <-
    (scores$combined_score - min) / (max - min)

  gene_names_to <-
    anno_df$alias[match(scores$to, anno_df$STRING_id)]
  gene_names_from <-
    anno_df$alias[match(scores$from, anno_df$STRING_id)]

  scores$to <- gene_names_to
  scores$from <- gene_names_from

  reverse_to <- scores$from
  reverse_from <- scores$to

  df <-
    data.frame(
      from = reverse_from,
      to = reverse_to,
      combined_score = scores$combined_score
    )
  scores <- dplyr::distinct(scores)
  df <- dplyr::distinct(df)

  scores <- rbind(scores, df)

  return(scores)
}
