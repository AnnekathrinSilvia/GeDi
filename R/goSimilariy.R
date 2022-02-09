#' Calculate the GO Similarity of all geneset combinations
#'
#' @param genesets a list of vectors of genes to score, one vector for each geneset
#' @param method the method to calculate the GO similarity see [GOSemSim] measure parameter
#' @param ontology the ontology of GO to use, see [GOSemSim] ont parameter
#' @param species the species of your data, indicated as org.XX.eg.db package from Bioconductor
#'
#' @return a [Matrix::Matrix()] with the pairwise GO similarity of each geneset pair
#' @export
#'
#' @examples
#' genesets <- list(c("PDHB", "VARS2", "IARS2"), c("IARS2", "PDHA2"))
#' go <- goSimilarity(genesets)
goSimilarity <- function(genesets, method = 'Wang', ontology = 'BP', species = 'org.Hs.eg.db'){
  l <- length(genesets)
  if(l == 0){
    return(-1)
  }
  go_sim <- Matrix::Matrix(0, l, l)
  go <- GOSemSim::godata(species, ont = ontology)

  for(i in 1:(l-1)){
    go1 <- genesets[[i]]
    for(j in (i+1):l){
      go2 <- genesets[[j]]
      sim <- GOSemSim::goSim(go1, go2, go, measure = method)
      go_sim[i, j] <- go_sim[j, i] <- sim
    }
  }
  return(go_sim)
}

#' Scaling computed distances with GO similarity scores
#'
#' @param scores a [Matrix::Matrix()] with the pairwise distances of each geneset pair
#' @param genesets a list of vectors of genes to score, one vector for each geneset
#' @param method the method to calculate the GO similarity see [GOSemSim] measure parameter
#' @param ontology the ontology of GO to use, see [GOSemSim] ont parameter
#' @param species the species of your data, indicated as org.XX.eg.db package from Bioconductor
#'
#' @return
#'
#' @examples
#' scores <- Matrix::Matrix(1, 1, 1)
#' genesets <- list(c("PDHB", "VARS2", "IARS2"), c("IARS2", "PDHA2"))
#' scaled <- scaleGO(scores, genesets)
scaleGO <- function(scores, genesets, method = 'Wang', ontology = 'BP', species = 'org.Hs.eg.db'){
  l <- nrow(scores)
  scaled <- Matrix::Matrix(0, l, l)
  scores_go <- goSimilarity(genesets, method, ontology, species)

  for(i in 1:(l-1)){
    for(j in (i+1):l){
      scaled[i, j] <- scaled[j, i] <- scores[i, j] * scores_go[i, j]
    }
  }
  return(scaled)
}
