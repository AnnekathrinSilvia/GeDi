#' Translates the species name into the NCBI ID
#'
#' @param species the species of the data to download
#'
#' @return the NCBI id of the species
#' @export
#' @importFrom taxize get_ids
#'
#' @examples
#' species <- "Homo sapiens"
#' getId(species)
getId <- function(species){
  id <- taxize::get_ids(species, db = "ncbi")
  return(id$ncbi[species])
}

#' Check if the downloaded PPI Matrix for a certain species needs to be updated
#'
#' @param bfc the location of the BioCFileCache
#' @param species the species of the data
#'
#' @return logical value if the PPI needs to be updated
#' @export
#' @importFrom BiocFileCache bfcquery bfcneedsupdate
#'
#' @examples
#' \dontrun{}
needsUpdate <- function(bfc, species){
  rid <- bfcquery(bfc, species)$rid
  return(bfcneedsupdate(bfc, rid))
}


#' Title
#'
#' @param species the species of the data
#' @param cachepath the path to save the cache
#' @param version the version of the StringDB database to download from
#'
#' @return the normalized PPI matrix
#' @export
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcadd bfcrpath
#'
#' @examples
#' \dontrun{}
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



#' List all PPIs saved in a Cache
#'
#' @param cachepath the path of the cache
#'
#' @return a list of PPIs
#' @export
#' @importFrom BiocFileCache BiocFileCache bfcinfo
#'
#' @examples
#' \dontrun{}
listPPI <- function(cachepath){
  bfc <- BiocFileCache(cachepath, ask = FALSE)
  bfc_df <- bfcinfo(bfc)

  return(bfc_df$rname)
}


#return as sparse matrix
#' Title
#'
#' @param ppi the path to a PPI matrix
#' @param info the path to the corresponding info file
#'
#' @return a normalized PPI dataframe
#' @export
#'
#' @examples
#' \dontrun{}
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
  colnames(ppi) <- columnnames

  return(ppi)
}


