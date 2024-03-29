#' Calculate similarity of GO terms
#'
#' Calculate the pairwise similarity of GO terms
#'
#' @param geneset_ids `list`, a `list` of GO identifiers to score
#' @param method character, the method to calculate the GO similarity.
#'               See [GOSemSim::goSim] measure parameter for possibilities.
#' @param ontology character, the ontology to use. See [GOSemSim::goSim]
#'                 `ont` parameter for possibilities.
#' @param species character, the species of your data. Indicated as
#'                org.XX.eg.db package from Bioconductor.
#' @param progress [shiny::Progress()] object, optional. To track the progress
#'                 of the function (e.g. in a Shiny app)
#' @param BPPARAM A BiocParallel `bpparam` object specifying how parallelization
#'                should be handled. Defaults to [BiocParallel::SerialParam()]
#'
#' @return A [Matrix::Matrix()] with the pairwise GO similarity of each
#'         geneset pair.
#' @export
#' @importFrom GOSemSim godata goSim
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom Matrix Matrix
#'
#' @examples
#'
#'
#' ## Mock example showing how the data should look like
#' go_ids <- c("GO:0002503", "GO:0045087", "GO:0019886",
#'             "GO:0002250", "GO:0001916", "GO:0019885")
#'
#' similarity <- goSimilarity(go_ids)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small, package = "GeDi")
#' go_ids <- macrophage_topGO_example_small$Genesets
#' \dontrun{
#' similarity <- goSimilarity(go_ids)
#' }
goSimilarity <- function(geneset_ids,
                         method = "Wang",
                         ontology = "BP",
                         species = "org.Hs.eg.db",
                         progress = NULL,
                         BPPARAM = BiocParallel::SerialParam()) {
  method <- match.arg(method, c("Resnik", "Lin", "Rel", "Jiang", "TCSS", "Wang"))

  if (method %in% c("Resnik", "Lin", "Rel", "Jiang"))
    useIC <- TRUE
  else
    useIC <- FALSE

  # Check if the species-specific org.XX.eg.db package is installed
  stopifnot("Species specific org.XX.eg.db
            is not installed" = system.file(package = species) != "")

   # Check if all geneset ids are GO identifiers
  go_ids <- all(vapply(geneset_ids, function(x) substr(x, 1, 2) == "GO",  logical(1)))
  stopifnot("Not all geneset ids are GO identifiers.
            This score only works on GO identifiers" = go_ids)

  # Determine the number of genesets
  l <- length(geneset_ids)
  if (l == 0) {
    return(-1)
  }

  # Initialize a matrix for GO similarity scores
  go_sim <- Matrix::Matrix(0, l, l)
  # Retrieve GO data for the specified species and ontology
  go <- godata(annoDb = species, ont = ontology, computeIC = useIC)

  results <- list()

  # Calculate the GO similarity for each pair of genesets
  for (g in seq_len((l - 1))) {
    a <- geneset_ids[[g]]
    if (!is.null(progress)) {
      progress$inc(1 / (l + 1), detail = paste("Scoring geneset number", g))
    }
    results[[g]] <- BiocParallel::bplapply((g + 1):l, function(i) {
      b <- geneset_ids[[i]]
      # Calculate GO similarity
      goSim(a, b, go, measure = method)
    }, BPPARAM = BPPARAM)
    go_sim[g, (g + 1):l] <- go_sim[(g + 1):l, g] <- unlist(results[[g]])
  }


  # Return the rounded GO similarity scores matrix
  return(round(go_sim, 2))
}

#' Scaling (distance) scores
#'
#' A method to scale a matrix of distance scores with the GO term similarity
#' of the associated genesets.
#'
#' @param scores a [Matrix::Matrix()], a matrix of (distance) scores for the
#'               identifiers in `geneset_ids`.
#' @param geneset_ids `list`, a `list` of GO identifiers to score
#' @param method character, the method to calculate the GO similarity.
#'               See [GOSemSim::goSim] measure parameter for possibilities.
#' @param ontology character, the ontology to use. See [GOSemSim::goSim]
#'                 ont parameter for possibilities.
#' @param species character, the species of your data. Indicated as
#'                org.XX.eg.db package from Bioconductor.
#' @param BPPARAM A BiocParallelParam object specifying how parallelization should
#'                be handled
#'
#' @return A [Matrix::Matrix()] of scaled values.
#' @export
#'
#' @examples
#' ## Mock example showing how the data should look like
#' go_ids <- c("GO:0002503", "GO:0045087", "GO:0019886",
#'             "GO:0002250", "GO:0001916", "GO:0019885")
#' set.seed(42)
#' scores <- Matrix::Matrix(stats::runif(36, min = 0, max = 1), 6, 6)
#' similarity <- scaleGO(scores,
#'                       go_ids)
#'
#' ## Example using the data available in the package
#' data(scores_macrophage_topGO_example_small, package = "GeDi")
#' data(macrophage_topGO_example_small, package = "GeDi")
#' go_ids <- macrophage_topGO_example_small$Genesets
#' \dontrun{
#' scores_scaled <- scaleGO(scores_macrophage_topGO_example_small,
#'                          go_ids)
#' }
scaleGO <- function(scores,
                    geneset_ids,
                    method = "Wang",
                    ontology = "BP",
                    species = "org.Hs.eg.db",
                    BPPARAM = BiocParallel::SerialParam()) {

  method <- match.arg(method, c("Resnik", "Lin", "Rel", "Jiang", "TCSS", "Wang"))

  # Determine the number of genesets
  l <- nrow(scores)

  # Ensure that the number of geneset_ids matches the number of genesets
  stopifnot(length(geneset_ids) == l)

  # Initialize a matrix for scaled scores
  scaled <- Matrix::Matrix(0, l, l)
  # Get GO similarity scores
  scores_go <- goSimilarity(geneset_ids, method, ontology, species, BPPARAM = BPPARAM)

  # Scale interaction scores with GO similarity scores
  for (i in seq_len((l - 1))) {
    for (j in (i + 1):l) {
      scaled[i, j] <- scaled[j, i] <- scores[i, j] * scores_go[i, j]
     }
  }

  # Return the scaled scores matrix
  return(scaled)
}
