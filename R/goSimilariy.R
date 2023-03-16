#' Calculate similarity of GO terms
#'
#' Calculate the pairwise similarity of GO terms
#'
#' @param geneset_ids `list`, a `list` of GO identifiers to score
#' @param method character, the method to calculate the GO similarity.
#'               See [GOSemSim::goSim] measure parameter for possibilities.
#' @param ontology character, the ontology to use. See [GOSemSim::goSim]
#'                 ont parameter for possibilities.
#' @param species character, the species of your data. Indicated as
#'                org.XX.eg.db package from Bioconductor.
#' @param progress [shiny::Progress()] object, optional. To track the progress
#'                 of the function (e.g. in a Shiny app)
#' @param n_cores numeric, number of cores to use for the function.
#'                Defaults to `NULL` in which case the function takes half of
#'                the available cores (see function .detectNumberCores(n_cores)).
#'
#' @return A [Matrix::Matrix()] with the pairwise GO similarity of each
#'         geneset pair.
#' @export
#' @import GOSemSim
#'
#' @examples
#' data(macrophage_topGO_example_small, package = "GeDi")
#' go_ids <- macrophage_topGO_example_small$Genesets
#' similarity <- goSimilarity(go_ids,
#'                            n_cores = 1)
goSimilarity <- function(geneset_ids,
                         method = "Wang",
                         ontology = "BP",
                         species = "org.Hs.eg.db",
                         progress = NULL,
                         n_cores = NULL) {
  stopifnot("Species specific org.XX.eg.db
            is not installed" = system.file(package = species) != "")
  go_ids <- all(sapply(geneset_ids, function(x) substr(x, 1, 2) == "GO"))
  stopifnot("Not all geneset ids are GO identifiers.
            This score only works on GO identifiers" = go_ids)
  l <- length(geneset_ids)
  if (l == 0) {
    return(-1)
  }
  go_sim <- Matrix::Matrix(0, l, l)
  go <- godata(species, ont = ontology)

  results <- list()

  # determine number of cores to use
  n_cores <- .getNumberCores(n_cores)

  # calculate the Jaccard distance for each pair of genesets
  for (g in 1:(l - 1)) {
    a <- geneset_ids[[g]]
    if (!is.null(progress)) {
      progress$inc(1 / (l + 1), detail = paste("Scoring geneset number", g))
    }
    results[[g]] <- parallel::mclapply((g + 1):l, function(i) {
      b <- geneset_ids[[i]]
      goSim(a, b, go, measure = method)
    }, mc.cores = n_cores)
    go_sim[g, (g + 1):l] <- go_sim[(g + 1):l, g] <- unlist(results[[g]])
  }
  #
  #   for (i in 1:(l - 1)) {
  #     go1 <- geneset_ids[[i]]
  #     if (!is.null(progress)) {
  #       progress$inc(1 / (l + 1), detail = paste("Scoring geneset number", i))
  #     }
  #     for (j in (i + 1):l) {
  #       go2 <- geneset_ids[[j]]
  #       sim <- goSim(go1, go2, go, measure = method)
  #       go_sim[i, j] <- go_sim[j, i] <- sim
  #     }
  #   }
  return(round(go_sim, 2))
}

#' Scaling (distance) scores
#'
#' A method to scale a matrix of distance scores with the GO term similarity
#' of the associated genesets.
#'
#' @param scores a [Matrix::Matrix()], a matrix of (distance) scores for the
#'               indentifiers in `geneset_ids`.
#' @param geneset_ids `list`, a `list` of GO identifiers to score
#' @param method character, the method to calculate the GO similarity.
#'               See [GOSemSim::goSim] measure parameter for possibilities.
#' @param ontology character, the ontology to use. See [GOSemSim::goSim]
#'                 ont parameter for possibilities.
#' @param species character, the species of your data. Indicated as
#'                org.XX.eg.db package from Bioconductor.
#' @param n_cores numeric, number of cores to use for the function.
#'                Defaults to `NULL` in which case the function takes half of
#'                the available cores (see function .detectNumberCores(n_cores)).
#'
#' @return Scaled Values
#' @export
#' @import GOSemSim
#'
#' @examples
#' data(scores_macrophage_topGO_example_small, package = "GeDi")
#' data(macrophage_topGO_example_small, package = "GeDi")
#' go_ids <- macrophage_topGO_example_small$Genesets
#' scores_scaled <- scaleGO(scores_macrophage_topGO_example_small,
#'                          go_ids,
#'                          n_cores = 1)
scaleGO <- function(scores,
                    geneset_ids,
                    method = "Wang",
                    ontology = "BP",
                    species = "org.Hs.eg.db",
                    n_cores = NULL) {
  l <- nrow(scores)
  stopifnot(length(geneset_ids) == l)
  scaled <- Matrix::Matrix(0, l, l)
  scores_go <- goSimilarity(geneset_ids, method, ontology, species, n_cores = n_cores)

  for (i in 1:(l - 1)) {
    for (j in (i + 1):l) {
      scaled[i, j] <- scaled[j, i] <- scores[i, j] * scores_go[i, j]
     }
  }
  return(scaled)
}
