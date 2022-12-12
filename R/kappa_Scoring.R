#' Calculate the Kappa distance
#'
#' Calculate the Kappa distance for two given genesets.
#'
#' @param a,b A vector of gene names whose interactions should
#'            be scored.
#' @param all_genes A vector of all unique genes in the data set from which `a`
#'                  and `b` are derived.
#'
#' @return The Kappa distance of the two genesets
#' @import dplyr
#' @export
#'
#' @examples
#' a <- c("PDHB", "VARS2")
#' b <- c("IARS2", "PDHA1")
#' all_genes <- c("PDHB", "VARS2", "IARS2", "PDHA1")
#' c <- calculateKappa(a, b, all_genes)
calculateKappa <- function(a, b, all_genes){
  n_genes <- length(all_genes)
  stopifnot(n_genes > 0)
  if(length(a) == 0 || length(b) == 0){
    return(1)
  }

  set_int <- length(intersect(a, b))
  set_none <- sum(!all_genes %in% a & !all_genes %in% b)

  only_a <- sum(all_genes %in% a & !all_genes %in% b)
  only_b <- sum(!all_genes %in% a & all_genes %in% b)

  total <- sum(only_a, only_b, set_int, set_none)

  O <- (set_int + set_none) / total
  E <- (set_int + only_a) * (set_int + only_b) + (only_b + set_none) * (only_a + set_none)
  E <- E / total^2

  kappa <- ((O-E) / (1-E))

  if(is.nan(kappa)){
    kappa <- 1
  }
  return(kappa)
}

#' Get Matrix of Kappa distances
#'
#' Calculate the Kappa distance of all combinations of genesets in a given data
#' set of genesets.
#'
#' @param genesets A `list` of genesets (each geneset is represented by a `list`
#'                 of the corresponding genes).
#' @param progress An optional [shiny::Progress()] object to track the progress
#'                 progress of the function in the app.
#'
#' @return A [Matrix::Matrix()] with the pairwise Kappa distance of each
#'         geneset pair.
#' @export
#'
#' @examples
#' genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
#' m <- getKappaMatrix(genesets)
getKappaMatrix <- function(genesets, progress = NULL){
  l <- length(genesets)
  if(l == 0){
    return(NULL)
  }
  k <- Matrix::Matrix(0, l, l)
  unique_genes <- unique(unlist(genesets))

  for(i in 1:(l - 1)){
    a <- genesets[[i]]
    for(j in (i+1):l){
      b <- genesets[[j]]
      k[i, j] <- k[j, i] <- calculateKappa(a, b, unique_genes)
    }
    if(!is.null(progress)){
      progress$inc(1/(l+1), detail = paste("Scoring geneset number", i))
    }
  }

  min <- min(k)
  max <- max(k)
  progress$inc(1/(l+1), detail = "Normalizing Kappa Matrix")
  for(i in 1:(l-1)){
    for(j in (i+1):l){
      k[i, j] <- k[j, i] <- 1 - ((k[i, j] - min)/(max - min))
    }
  }
  return(k)
}
