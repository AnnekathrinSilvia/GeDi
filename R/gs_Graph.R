#' Title
#'
#' @param distanceMatrix a Distance Matrix containing for each pair the distance
#' @param cutOff a cutOff value classifying a threshold for similarity
#'
#' @return
#' @import Matrix
#' @export
#'
#' @examples
getAdjacencyMatrix <- function(distanceMatrix, cutOff){
  l <- nrow(distanceMatrix)
  adjMat <- Matrix::Matrix(0, l, l)

  for(i in 1:l){
    edge <- which(distanceMatrix[i, ] <= cutOff)

    adjMat[i, edge] <- 1
  }
  rownames(adjMat) <- rownames(distanceMatrix)
  colnames(adjMat) <- colnames(distanceMatrix)
  return(adjMat)
}

#' Title
#'
#' @param adjMatrix
#'
#' @return
#' @import igraph
#' @export
#'
#' @examples
buildGraph <- function(adjMatrix){
  g <- igraph::graph_from_adjacency_matrix(
    adjMatrix,
    mode = "undirected",
    add.colnames = NULL,
    add.rownames = NA,
    diag = F
  )

  V(g)$color <- "#0092AC"

  return(g)
}

#' Title
#'
#' @param cluster
#' @param geneset_names
#'
#' @return
#' @export
#'
#' @examples
getClusterAdjacencyMatrix <- function(cluster, geneset_names){
  l <- length(geneset_names)
  adj <- Matrix::Matrix(0, l, l)

  for(i in 1:length(cluster)){
    subcluster <- cluster[[i]]
    adj[subcluster, subcluster] <- 1
  }

  rownames(adj) <-  colnames(adj) <- geneset_names
  diag(adj) <- 0
  return(adj)
}

#' Title
#'
#' @param cluster
#' @param geneset_names
#'
#' @return
#' @export
#'
#' @examples
getBipartiteGraph <- function(cluster, geneset_names){
  edgelist <- c()
  type <- c()
  n_cluster <- length(cluster)
  node_number <- n_cluster + 1
  df_node_mapping <- data.frame(matrix(NA, nrow = length(geneset_names), ncol = 1))

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

  V(graph)$nodeType <- NA
  V(graph)$nodeType[cluster_id] <- "Cluster"
  V(graph)$nodeType[geneset_id] <- "Geneset"

  V(graph)$shape <- c("box", "ellipse")[factor(V(graph)$nodeType, levels = c("Cluster", "Geneset"))]

  V(graph)$color[cluster_id] <- "gold"
  V(graph)$color[geneset_id] <- "#0092AC"
  E(graph)$color <- "black"

  return(graph)
}
