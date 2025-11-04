#' Calculate similarity of GO terms
#'
#' Calculate the pairwise similarity of GO terms
#'
#' @param geneset_ids `list`, a `list` of GO identifiers to score
#' @param method character, the method to calculate the GO distance.
#'               Possible options are "Resnik", "Lin", "Rel",
#'               "Jiang",  "Wang".
#' @param ontology character, the ontology to use. Possible options are "BP",
#'                 "MF" and "CC".
#' @param species character, the species of your data. Indicated as
#'                org.XX.eg.db package from Bioconductor.
#'
#' @return A [Matrix::Matrix()] with the pairwise GO distance of each
#'         geneset pair.
#' @export
#' @importFrom simona term_sim create_ontology_DAG_from_GO_db dag_has_terms
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
  dag <- create_ontology_DAG_from_GO_db(ontology, org_db = species)

  if (method == "Resnik") {
    method <- "Sim_Resnik_1999"
  } else if (method == "Lin") {
    method <- "Sim_Lin_1998"
  } else if (method == "Rel") {
    method <- "Sim_Relevance_2006"
  } else if (method == "Jiang") {
    method <- "Sim_Jiang_1997"
  } else if (method == "Wang") {
    method <- "Sim_Wang_2007"
  }

  if (method == "Sim_Resnik_1999") {
    go_sim <- term_sim(dag = dag,
                       terms = geneset_ids,
                       method = method,
                       control = list(norm_method = "Nunif"))
  } else if (method == "Sim_Jiang_1997") {
    go_sim <- term_sim(dag = dag,
                       terms = geneset_ids,
                       method = method,
                       control = list(norm_method = "max"))
  } else {
    go_sim <- term_sim(dag = dag,
                       terms = geneset_ids,
                       method = method)
  }

  terms_in_dag <- dag_has_terms(dag, geneset_ids)
  ids <- geneset_ids[terms_in_dag]
  rownames(go_sim) <- colnames(go_sim) <- ids

  go_dist <- as.matrix(1 - go_sim)
  diag(go_dist) <- 0
  # Return the rounded GO distance scores matrix
  return(round(go_dist, 2))
}
