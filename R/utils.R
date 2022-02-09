#Splits the genesets into lists of genes, one for each gene set
#' Split the genesets into a list of vectors of genes, one for each geneset
#'
#' @param genesets a list of the genesets
#'
#' @return a list of vectors of genes, on for each geneset
#' @export
#'
#' @examples
#' genesets <- list(c("PDHB", "VARS2", "IARS2"), c("IARS2", "PDHA2"))
#' genes <- getGenes(genesets)
getGenes <- function(genesets){
  if(length(genesets) == 0){
    return(list())
  }
  genesets <- lapply(1:nrow(genesets), function(i) {
    toupper(strsplit(genesets[i, 2], " ")[[1]])
  })

  return(genesets)
}


#' Assign each gene in the list to its index in the PPI
#'
#' @param genes a list of vectors of genes
#' @param ppi a Protein-Protein interaction matrix
#'
#' @return a list of vectors of indices
#' @export
#'
#' @examples
#' genes <- list(c("PDHB", "VARS2", "IARS2"), c("IARS2", "PDHA2"))
#' ppi <- Matrix::matrix(0.5, 4, 4)
#' genes_indexed <- getGenesIndexed(genes, ppi)
getGenesIndexed <- function(genes, ppi){
  if(length(genes) == 0){
    return(list())
  }
  genes <- sapply(1:length(genes), function(i) {
    unname(unlist(sapply(genes[[i]],
                         function(j) {
                           a <- which(j == rownames(ppi))
                           if (length(a)) {
                             a
                           }
                         }), use.names = FALSE))
  })
  return(genes)
}
