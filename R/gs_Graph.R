#' Construct an adjacency matrix
#'
#' Construct an adjacency matrix from the (distance) scores and a given
#' threshold.
#'
#' @param distanceMatrix A [Matrix::Matrix()] containing (distance) scores
#'                       between 0 and 1.
#' @param cutOff Numeric value, indicating for which pair of entries in the
#'               `distanceMatrix` a 1 should be inserted in the adjacency matrix.
#'               A 1 is inserted when for each entry in the matrix that is
#'               smaller or equal to the `cutOff` value.
#'
#' @return A [Matrix::Matrix()] of adjacency status
#' @importFrom Matrix Matrix
#' @export
#'
#' @examples
#' m <- Matrix::Matrix(runif(1000, 0, 1), 100, 100)
#' threshold <- 0.3
#' adj <- getAdjacencyMatrix(m, threshold)
getAdjacencyMatrix <- function(distanceMatrix, cutOff){
  if(is.null(distanceMatrix)){
    return(NULL)
  }
  l <- nrow(distanceMatrix)
  adjMat <- Matrix::Matrix(0, l, l)

  for(i in 1:l){
    edge <- which(distanceMatrix[i, ] <= cutOff)
    if(length(edge) > 0){
      adjMat[i, edge] <- 1
    }
  }
  rownames(adjMat) <- rownames(distanceMatrix)
  colnames(adjMat) <- colnames(distanceMatrix)
  return(adjMat)
}

#' Construct a graph
#'
#' Construct a graph from a given adjacency matrix
#'
#' @param adjMatrix A [Matrix::Matrix()] indicating for which pair of nodes an
#'                  edge should be added; 1 indicating an edge, 0 indicating no
#'                  edge.
#'
#' @return An `igraph` object to be further manipulated or processed/plotted
#'         (e.g. via [igraph::plot.igraph()] or
#'         [visNetwork::visIgraph()][visNetwork::visNetwork-igraph])
#' @importFrom igraph graph_from_adjacency_matrix V
#' @export
#'
#' @examples
#'
#' adj <- Matrix::Matrix(0, 100, 100)
#' adj[c(80:100), c(80:100)] <- 1
#' graph <- buildGraph(adj)
buildGraph <- function(adjMatrix){
  g <- igraph::graph_from_adjacency_matrix(
    adjMatrix,
    mode = "undirected",
    add.colnames = NULL,
    add.rownames = NA,
    diag = F
  )

  igraph::V(g)$color <- "#0092AC"

  return(g)
}

#' Construct an adjacency matrix
#'
#' Construct an adjacency matrix from a `list` of cluster..
#'
#' @param cluster A `list` of clusters, where each cluster member is indicated
#'                by a numeric value
#' @param geneset_names A vector of geneset names
#'
#' @return A [Matrix::Matrix()] of adjacency status
#' @importFrom Matrix Matrix
#' @export
#'
#' @examples
#' cluster <- list(c(1:5), c(6:9))
#' geneset_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
#' adj <- getClusterAdjacencyMatrix(cluster, geneset_names)
getClusterAdjacencyMatrix <- function(cluster, geneset_names){
  l <- length(geneset_names)
  adj <- Matrix::Matrix(0, l, l)
  if(length(cluster) == 0){
    rownames(adj) <-  colnames(adj) <- geneset_names
    diag(adj) <- 0
    return(adj)
  }

  for(i in 1:length(cluster)){
    subcluster <- cluster[[i]]
    adj[subcluster, subcluster] <- 1
  }

  rownames(adj) <- colnames(adj) <- geneset_names
  diag(adj) <- 0
  return(adj)
}

#' Construct a bipartite graph
#'
#' Construct a bipartite graph from cluster information, mapping the cluster
#' to its members
#'
#' @param cluster A `list` clusters, where each cluster member is indicated
#'                by a numeric value
#' @param geneset_names A vector of geneset names
#' @param genes A `list` of `list` of genes which belong to the genesets in geneset_names
#'
#' @return An `igraph` object to be further manipulated or processed/plotted
#'         (e.g. via [igraph::plot.igraph()] or
#'         [visNetwork::visIgraph()][visNetwork::visNetwork-igraph])
#' @export
#' @importFrom igraph make_bipartite_graph set_vertex_attr V E
#' @importFrom stats na.omit
#'
#' @examples
#' cluster <- list(c(1:5), c(6:9))
#' geneset_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
#' genes <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
getBipartiteGraph <- function(cluster, geneset_names, genes){
  stopifnot(length(cluster) > 0)
  stopifnot(length(geneset_names) > 0)
  stopifnot(length(genes) > 0)

  edgelist <- c()
  type <- c()
  n_cluster <- length(cluster)
  node_number <- n_cluster + 1
  df_node_mapping <- data.frame(matrix(NA, nrow = length(geneset_names), ncol = 1))
  colnames(df_node_mapping) <- "Node_number"

  node_labels <- c()

  for(i in 1:n_cluster){
    node_labels <- c(node_labels, paste0("Cluster ", i))
  }

  for(i in 1:n_cluster){
    subcluster <- cluster[[i]]
    for(j in subcluster){
      if(!is.na(df_node_mapping[j, ])){
        n <- df_node_mapping[j, ]
        edgelist <- c(edgelist, i, n)
      }else{
        edgelist <- c(edgelist, i, node_number)
        df_node_mapping[j, ] <- node_number
        node_number <- node_number + 1
        node_labels <- c(node_labels, geneset_names[[j]])
      }
    }
  }

  type <- c(rep(0, n_cluster))
  type <- c(type, rep(1, node_number-n_cluster-1))

  graph <- igraph::make_bipartite_graph(type, edgelist, directed = TRUE)
  graph <- set_vertex_attr(graph, "name", value = node_labels)

  cluster_id <- which(names(V(graph)) %in% node_labels[1:n_cluster])
  geneset_id <- which(!(names(V(graph)) %in% node_labels[1:n_cluster]))

  igraph::V(graph)$nodeType <- NA
  igraph::V(graph)$nodeType[cluster_id] <- "Cluster"
  igraph::V(graph)$nodeType[geneset_id] <- "Geneset"

  igraph::V(graph)$shape <- c("box", "ellipse")[factor(V(graph)$nodeType, levels = c("Cluster", "Geneset"))]

  igraph::V(graph)$color[cluster_id] <- "gold"
  igraph::V(graph)$color[geneset_id] <- "#0092AC"
  igraph::E(graph)$color <- "black"


  igraph::V(graph)$title <- NA

  for(i in cluster_id){
    igraph::V(graph)$title[i] <- paste0(
      "<h4>", igraph::V(graph)$name[i], "</h4><br>",
      "Members = ", paste(geneset_names[(cluster[[i]])], collapse = ";")
    )
  }

  for(j in geneset_id){
    igraph::V(graph)$title[j] <- paste0(
      "<h4>", igraph::V(graph)$name[j], "</h4><br>",
      "Members = ", paste(genes[as.integer(stats::na.omit(rownames(df_node_mapping)[df_node_mapping$Node_number == j]))], collapse = " ")
    )
  }

  return(graph)
}
