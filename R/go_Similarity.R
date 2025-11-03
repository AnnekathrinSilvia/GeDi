#' Calculate similarity of GO terms
#'
#' Calculate the pairwise similarity of GO terms
#'
#' @param geneset_ids `list`, a `list` of GO identifiers to score
#' @param method character, the method to calculate the GO distance.
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
#' @return A [Matrix::Matrix()] with the pairwise GO distance of each
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
#' similarity <- goDistance(go_ids)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small, package = "GeDi")
#' go_ids <- macrophage_topGO_example_small$Genesets
#' \dontrun{
#' similarity <- goDistance(go_ids)
#' }
goDistance <- function(geneset_ids,
                         method = "Wang",
                         ontology = "BP",
                         species = "org.Hs.eg.db",
                         progress = NULL,
                         BPPARAM = BiocParallel::SerialParam()) {
  method <- match.arg(method, c("Resnik", "Lin", "Rel",
                                "Jiang", "TCSS", "Wang"))
  if (method %in% c("Resnik", "Lin", "Rel", "Jiang"))
    useIC <- TRUE
  else
    useIC <- FALSE

  # Check if the species-specific org.XX.eg.db package is installed
  stopifnot("Species specific org.XX.eg.db
            is not installed" = system.file(package = species) != "")
  # Check if all geneset ids are GO identifiers
  go_ids <- all(vapply(geneset_ids, function(x) substr(x, 1, 2) == "GO",
                       logical(1)))
  stopifnot("Not all geneset ids are GO identifiers.
            This score only works on GO identifiers" = go_ids)
  # Determine the number of genesets
  l <- length(geneset_ids)
  if (l == 0) {
    return(-1)
  }

  # Initialize a matrix for GO distance scores
  go_sim <- Matrix::Matrix(0, l, l)
  # Retrieve GO data for the specified species and ontology
  go <- godata(annoDb = species, ont = ontology, computeIC = useIC)

  results <- list()
  # Calculate the GO distance for each pair of genesets
  for (g in seq_len((l - 1))) {
    a <- geneset_ids[[g]]
    if (!is.null(progress)) {
      progress$inc(1 / (l + 1), detail = paste("Scoring geneset number", g))
    }
    results[[g]] <- BiocParallel::bplapply((g + 1):l, function(i) {
      b <- geneset_ids[[i]]
      # Calculate GO distance
      goSim(a, b, go, measure = method)
    }, BPPARAM = BPPARAM)
    go_sim[g, (g + 1):l] <- go_sim[(g + 1):l, g] <- unlist(results[[g]])
  }
  
  go_sim[is.na(go_sim)] <- 0
  
  # Next, we have to normalize some of the similarities to the [0, 1]
  # interval and afterwards transform the similarity to a distance by 
  # calculating 1 - Similarity
  if(method %in% c("Resnik", "Jiang")){
    min <- min(go_sim)
    max <- max(go_sim)
    
    if (!is.null(progress)) {
      progress$inc(1 / (l + 1), detail = "Normalizing Similarity Matrix")
    }
    
    # Update the matrix with normalized values
    results <- list()
    for (j in seq_len((l - 1))) {
      results[[j]] <- bplapply((j + 1):l, function(i) {
        return(1 - ((go_sim[j, i] - min) / (max - min)))
      }, BPPARAM = BPPARAM)
      go_sim[j, (j + 1):l] <- go_sim[(j + 1):l, j] <- unlist(results[[j]])
    }
  }
  
  go_dist <- 1 - go_sim
  diag(go_dist) <- 0
  # Return the rounded GO distance scores matrix
  return(round(go_dist, 2))
}



#' Calculate similarity of GO terms
#'
#' Calculate the pairwise similarity of GO terms
#'
#' @param geneset_ids `list`, a `list` of GO identifiers to score
#' @param method character, the method to calculate the GO distance.
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
#' @return A [Matrix::Matrix()] with the pairwise GO distance of each
#'         geneset pair.
#' @export
#' @importFrom simona term_sim create_ontology_DAG_from_GO_db
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
#' similarity <- goDistance_(go_ids)
#'
#' ## Example using the data available in the package
#' data(macrophage_topGO_example_small, package = "GeDi")
#' go_ids <- macrophage_topGO_example_small$Genesets
#' \dontrun{
#' similarity_revision <- goDistance_REVISION(go_ids)
#' }
goDistance_REVISION <- function(geneset_ids,
                       method = "Wang",
                       ontology = "BP",
                       species = "org.Hs.eg.db") {
  method <- match.arg(method, c("Resnik", "Lin", "Rel",
                                "Jiang",  "Wang"))

  # Check if the species-specific org.XX.eg.db package is installed
  stopifnot("Species specific org.XX.eg.db
            is not installed" = system.file(package = species) != "")
  # Check if all geneset ids are GO identifiers
  go_ids <- all(vapply(geneset_ids, function(x) substr(x, 1, 2) == "GO",
                       logical(1)))
  stopifnot("Not all geneset ids are GO identifiers.
            This score only works on GO identifiers" = go_ids)
  # Determine the number of genesets
  l <- length(geneset_ids)
  if (l == 0) {
    return(-1)
  }
  
  # Initialize a matrix for GO distance scores
  go_sim <- Matrix::Matrix(0, l, l)
  
  # Create the dag
   dag = create_ontology_DAG_from_GO_db(ontology, org_db = species)

  if(method == "Resnik"){
    method = "Sim_Resnik_1999"
  }else if (method == "Lin"){
      method = "Sim_Lin_1998"
  }else if(method == "Rel"){
      method = "Sim_Relevance_2006"
  }else if(method == "Jiang"){
    method = "Sim_Jiang_1997"
  }else if(method == "Wang"){
    method = "Sim_Wang_2007"
  }

  if(method == "Sim_Resnik_1999"){
    go_sim <- term_sim(dag = dag,
                       terms = geneset_ids,
                       method = method,
                       control = list(norm_method = "Nunif"))
  } else if(method == "Sim_Jiang_1997"){
    go_sim <- term_sim(dag = dag,
                       terms = geneset_ids,
                       method = method,
                       control = list(norm_method = "max"))
  }else{
    go_sim <- term_sim(dag = dag,
                       terms = geneset_ids,
                       method = method)
  }

  go_dist <- as.matrix(1 - go_sim)
  diag(go_dist) <- 0
  # Return the rounded GO distance scores matrix
  return(round(go_dist, 2))
}

### IMPLEM: the corresponding benchmark steps -----
# library("fastmatch")
# library("proxyC")
# 
# microbenchmark::microbenchmark(
#   goDistance(go_ids),
#   goDistance_REVISION(go_ids),
#   times = 100
# )
# 
# bench::mark(
#   goDistance(go_ids),
#   goDistance_REVISION(go_ids),
#   iterations = 100,
#   memory = FALSE,
#   check = FALSE
# )

