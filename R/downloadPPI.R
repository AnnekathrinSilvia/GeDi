#' Title
#'
#' @param species
#'
#' @return
#' @export
#'
#' @examples
getId <- function(species){
  id <- taxize::get_ids(species, db = "ncbi")
  return(id$ncbi[species])
}

#' Title
#'
#' @param bfc
#' @param species
#'
#' @return
#' @export
#' @importFrom BiocFileCache bfcquery bfcneedsupdate
#'
#' @examples
needsUpdate <- function(bfc, species){
  rid <- bfcquery(bfc, species)$rid
  return(bfcneedsupdate(bfc, rid))
}


#' Title
#'
#' @param species
#' @param cachepath
#' @param version
#'
#' @return
#' @export
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcadd bfcrpath
#'
#' @examples
downloadPPI <- function(species, cachepath = "~/GSD", version = "11.5"){
  bfc <- BiocFileCache(cachepath, ask = FALSE)
  name_ppi <- paste("PPI_",
                    species,
                    sep = "")

  name_info <- paste("PPI_Info_",
                     species,
                     sep = "")

  if(nrow(bfcquery(bfc, name_ppi)) == 0 || needsUpdate(bfc, name_ppi)){
    species_id <- getId(species)
    url_ppi <- paste("https://stringdb-static.org/download/protein.links.v",
                     version,
                     "/",
                     species_id,
                     ".protein.links.v",
                     version,
                     ".txt.gz",
                     sep = "")

    url_info <- paste("https://stringdb-static.org/download/protein.info.v",
                      version,
                      "/",
                      species_id,
                      ".protein.info.v",
                      version,
                      ".txt.gz",
                      sep = "")

    ppi <- bfcadd(bfc, rname = name_ppi, fpath = url_ppi)
    info <- bfcadd(bfc, rname = name_info, fpath = url_info)

    ppi <- normalizePPI(ppi, info)
    return(ppi)
  }else{
  ppi <- bfcrpath(bfc, name_ppi)
  info <- bfcrpath(bfc, name_info)

  ppi <- normalizePPI(ppi, info)

  return(ppi)
  }
}



# libPaths to get to where libraries are installed, then put them in the gsd folder under PPIs
#-> not good maybe users have a strange way of installing every package somewhere else
# -> rather use "C:\\Users\\anneludt\\AppData\\Local/R/cache/R/AnnotationHub
#    (with GSD/package name instead of AnnotationHub)



#' Title
#'
#' @param cachepath
#'
#' @return
#' @export
#'
#' @examples
listPPI <- function(cachepath){
  bfc <- BiocFileCache(cachepath, ask = FALSE)
  bfc_df <- bfcinfo(bfc)

  return(bfc_df$rname)
}


#return as sparse matrix
#' Title
#'
#' @param ppi
#' @param info
#'
#' @return
#' @export
#'
#' @examples
normalizePPI <- function(ppi, info){
  ppi <- as.data.frame(read.delim(ppi, sep = " ", header = TRUE))
  info <- as.data.frame(read.delim(info, sep = "\t", fill = FALSE, quote = ""))

  max <- max(ppi$combined_score)
  min <- min(ppi$combined_score)
  rownames <- info$preferred_name[match(ppi$protein1, info$X.string_protein_id)]
  columnnames <- info$preferred_name[match(ppi$protein2, info$X.string_protein_id)]

  scores <- list()
  new_scores <- lapply(ppi$combined_score, function(x) ((x - min) / (max - min)))


  ppi <- as.data.frame(new_scores)
  rownames(ppi) <- rownames
  columnnames(ppi) <- columnnames

  return(ppi)
}


