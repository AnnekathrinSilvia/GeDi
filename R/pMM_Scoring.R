#' Sum up Protein-Protein interactions between two genesets
#'
#' The function sums up the strength of the Protein-Protein interactions of
#' all genes in the given genesets.
#'
#' @param a,b A vector of gene names whose interactions should
#'            be scored.
#' @param ppi A `data.frame` object which contains the Protein-Protein
#'            interactions.
#'            The object has three columns, `from` and `to` which
#'            specify the gene names of the interacting proteins in no
#'            particular order (symmetric interaction) and a column
#'            `combined_score` which is a numerical value of the strength of
#'            the interaction.
#'
#'
#' @return The sum of the Protein-Protein interactions.
#' @export
#' @import dplyr
#' @examples
#' a <- c("PDHB", "VARS2")
#' b <- c("IARS2", "PDHA1")
#'
#' ppi <- data.frame(from = c("PDHB", "VARS2"),
#'                   to = c("IARS2", "PDHA1"),
#'                   combined_score = c(0.5, 0.2))
#'
#' score <- sumInteraction(a, b, ppi)
sumInteraction <- function(a, b, ppi) {
  if (length(a) ==  0 | length(b) == 0) {
    return(0)
  } else{
    interactions <- which(ppi$from %in% a & ppi$to %in% b)
    return(sum(ppi[interactions,]$combined_score))
  }
}

#' Calculate interaction score for two genesets
#'
#' The function calculates the interaction score for two given genesets.
#'
#' @param a,b A vector of gene names whose interactions should
#'            be used.
#' @param ppi A `data.frame` object which contains the Protein-Protein
#'            interactions.
#'            The object has three columns, `from` and `to` which
#'            specify the gene names of the interacting proteins in no
#'            particular order (symmetric interaction) and a column
#'            `combined_score` which is a numerical value of the strength of
#'            the interaction.
#' @param maxInteract The maximum value of the `combined_score` column of
#'                    the `ppi`.
#'
#' @return The interaction score of the two given genesets.
#' @export
#'
#' @examples
#' a <- c("PDHB", "VARS2")
#' b <- c("IARS2", "PDHA1")
#'
#' ppi <- data.frame(from = c("PDHB", "VARS2"),
#'                   to = c("IARS2", "PDHA1"),
#'                   combined_score = c(0.5, 0.2))
#' maxInteract <- max(ppi$combined_score)
#'
#' interaction <- getInteractionScore(a, b, ppi, maxInteract)
getInteractionScore <- function(a, b, ppi, maxInteract) {
  if(length(a) == 0 || length(b) == 0){
    return(0)
  }
  onlya <- setdiff(a, b)
  onlyb <- setdiff(b, a)
  int <- intersect(a, b)

  w <- min(length(a), length(b)) / (length(a) + length(b))
  intlength <- length(int)
  onlyblength <- length(onlyb)

  sumInt <- sumInteraction(onlya, int, ppi)
  sumOnlyb <- sumInteraction(onlya, onlyb, ppi)

  nom <- (w * sumInt) + sumOnlyb
  denom <- maxInteract * (w * intlength + onlyblength)

  return(nom / denom)
}

#' Calculate local pMM distance
#'
#' Calculate the local pMM distance of two genesets.
#'
#' @param a,b A vector of gene names whose interactions should
#'            be scored.
#' @param ppi A `data.frame` object which contains the Protein-Protein
#'            interactions.
#'            The object has three columns, `from` and `to` which
#'            specify the gene names of the interacting proteins in no
#'            particular order (symmetric interaction) and a column
#'            `combined_score` which is a numerical value of the strength of
#'            the interaction.
#' @param maxInteract The maximum value of the `combined_score` column of
#'                    the `ppi`.
#' @param alpha A scaling factor (between 0 and 1) which indicates how much the
#'              Protein-Protein interactions are weighted into the final score.
#'              Defaults to 1.
#'
#' @return The pMM distance for two genesets.
#' @export
#'
#' @examples
#' a <- c("PDHB", "VARS2")
#' b <- c("IARS2", "PDHA1")
#'
#' ppi <- data.frame(from = c("PDHB", "VARS2"),
#'                   to = c("IARS2", "PDHA1"),
#'                   combined_score = c(0.5, 0.2))
#' maxInteract <- max(ppi$combined_score)
#'
#' pMM_score <- pMMlocal(a, b, ppi, maxInteract)
pMMlocal <- function(a, b, ppi, maxInteract, alpha = 1) {
  z <- min(length(a), length(b))
  if(z == 0){
    return(1)
  }

  factor1 <- (length(intersect(a, b))) / z
  factor2 <- (alpha / z) * getInteractionScore(a, b, ppi, maxInteract)

  return(min(factor1 + factor2, 1))
}


#' Calculate the pMM distance
#'
#' Calculate the pMM distance of all combinations of genesets in a given data
#' set of genesets.
#'
#' @param genes A `list` of genesets (each geneset is represented by a vector
#'              of the corresponding genes).
#' @param ppi A `data.frame` object which contains the Protein-Protein
#'            interactions.
#'            The object has three columns, `from` and `to` which
#'            specify the gene names of the interacting proteins in no
#'            particular order (symmetric interaction) and a column
#'            `combined_score` which is a numerical value of the strength of
#'            the interaction.
#' @param alpha A scaling factor (between 0 and 1) which indicates how much the
#'              Protein-Protein interactions are weighted into the final score.
#'              Defaults to 1.
#' @param progress An optional [shiny::Progress()] object to track the progress
#'                 progress of the function in the app.
#'
#' @return A [Matrix::Matrix()] with the pairwise pMM distance of each
#'         geneset pair.
#' @export
#'
#' @examples
#' genes <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
#'
#' ppi <- data.frame(from = c("PDHB", "VARS2"),
#'                   to = c("IARS2", "PDHA1"),
#'                   combined_score = c(0.5, 0.2))
#'
#' pMM <- getpMMMatrix(genes, ppi)
getpMMMatrix <- function(genes, ppi, alpha = 1, progress = NULL){
  l <- length(genes)
  if(l == 0){
    return(NULL)
  }
  stopifnot(!is.null(ppi))
  scores <- Matrix::Matrix(0, l, l)
  maxInteract <- max(ppi$combined_score)
  for (i in 1:(l - 1)) {
    a <- genes[[i]]
    for (j in (i + 1):l) {
      b <- genes[[j]]
      pmm <- min(
        pMMlocal(a, b, ppi, maxInteract),
        pMMlocal(b, a, ppi, maxInteract)
      )
      scores[i, j] <- scores[j, i] <- pmm
    }
    if(!is.null(progress)){
      progress$inc(1/l, detail = paste("Scoring geneset number", i))
    }
  }
  return(scores)
}
