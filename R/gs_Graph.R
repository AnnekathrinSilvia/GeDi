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
#' @import igraph
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

  gs_names <- rownames(adjMatrix)
  ids <- which(names(V(g)) %in% gs_names)

  V(g)$title[ids] <- paste0(
    "<h4>",
    sprintf('<a href="http://amigo.geneontology.org/amigo/term/%s" target="_blank">%s</a>', gs_names[ids], gs_names[ids]), "</h4><br>",
    V(g)$name[ids], "<br><br>")

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

#' Title
#'
#' @param cluster
#' @param geneset_df
#' @param gs_names
#' @param color_by
#'
#' @return
#' @import igraph
#' @importFrom GeneTonic map2color
#' @importFrom grDevices colorRampPalette
#' @export
#'
#' @examples
buildClusterGraph <- function(cluster,
                              geneset_df,
                              gs_names,
                              color_by = NULL){

  adj <- getClusterAdjacencyMatrix(cluster,
                                   gs_names)
  g <- buildGraph(adj)

  if(!is.null(color_by)){
    if (!color_by %in% colnames(geneset_df)) {
      stop(
        "Your data does not contain the ",
        color_by,
        " column.\n",
        "Please select another column to use for the color."
      )
    }
  }

  ids <- which(names(V(g)) %in% gs_names)

  if(!is.null(color_by)){
    col_var <- geneset_df[ids, color_by]
    # the palette changes if it is z_score VS pvalue
    if (all(col_var <= 1) & all(col_var > 0)) { # likely p-values...
      col_var <- -log10(col_var)
      # V(g)$color <- colVar
      mypal <- (scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "YlOrRd", 9))(50), 0.8
      ))
      mypal_hover <- (scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "YlOrRd", 9))(50), 0.5
      ))
      mypal_select <- (scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "YlOrRd", 9))(50), 1
      ))

      V(g)$color.background <- map2color(col_var, mypal, symmetric = FALSE,
                                         limits = range(na.omit(col_var)))
      V(g)$color.highlight <- map2color(col_var, mypal_select, symmetric = FALSE,
                                        limits = range(na.omit(col_var)))
      V(g)$color.hover <- map2color(col_var, mypal_hover, symmetric = FALSE,
                                    limits = range(na.omit(col_var)))

      V(g)$color.background[is.na(V(g)$color.background)] <- "lightgrey"
      V(g)$color.highlight[is.na(V(g)$color.highlight)] <- "lightgrey"
      V(g)$color.hover[is.na(V(g)$color.hover)] <- "lightgrey"
    } else {
      # e.g. using z_score or aggregated value
      if (prod(range(na.omit(col_var))) >= 0) {
        # gradient palette
        mypal <- (scales::alpha(
          colorRampPalette(RColorBrewer::brewer.pal(name = "Oranges", 9))(50), 0.8
        ))
        mypal_hover <- (scales::alpha(
          colorRampPalette(RColorBrewer::brewer.pal(name = "Oranges", 9))(50), 0.5
        ))
        mypal_select <- (scales::alpha(
          colorRampPalette(RColorBrewer::brewer.pal(name = "Oranges", 9))(50), 1
        ))

        V(g)$color.background <- map2color(col_var, mypal, symmetric = FALSE,
                                           limits = range(na.omit(col_var)))
        V(g)$color.highlight <- map2color(col_var, mypal_select, symmetric = FALSE,
                                          limits = range(na.omit(col_var)))
        V(g)$color.hover <- map2color(col_var, mypal_hover, symmetric = FALSE,
                                      limits = range(na.omit(col_var)))
        V(g)$color.background[is.na(V(g)$color.background)] <- "lightgrey"
        V(g)$color.highlight[is.na(V(g)$color.highlight)] <- "lightgrey"
        V(g)$color.hover[is.na(V(g)$color.hover)] <- "lightgrey"

      } else {
        # divergent palette to be used
        mypal <- rev(scales::alpha(
          colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.8
        ))
        mypal_hover <- rev(scales::alpha(
          colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.5
        ))
        mypal_select <- rev(scales::alpha(
          colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 1
        ))

        V(g)$color.background <- map2color(col_var, mypal, symmetric = TRUE,
                                           limits = range(na.omit(col_var)))
        V(g)$color.highlight <- map2color(col_var, mypal_select, symmetric = TRUE,
                                          limits = range(na.omit(col_var)))
        V(g)$color.hover <- map2color(col_var, mypal_hover, symmetric = TRUE,
                                      limits = range(na.omit(col_var)))

        V(g)$color.background[is.na(V(g)$color.background)] <- "lightgrey"
        V(g)$color.highlight[is.na(V(g)$color.highlight)] <- "lightgrey"
        V(g)$color.hover[is.na(V(g)$color.hover)] <- "lightgrey"
      }
    }
  }

  no_cluster <- V(g)[degree(g) == 0]
  g <- delete_vertices(g, no_cluster)

  transposed_df <- as.data.frame(t(geneset_df))


  title <- list()
  names_rows <- rownames(transposed_df)

  for(i in 1:ncol(transposed_df)){
    node_title <-"<!DOCTYPE html> <html> <head> <style>
      table {font-family: arial, sans-serif; font-size: 10px; border-collapse: collapse;width: 100%;} td,
      th { border: 1px solid #dddddd; text-align: center; padding: 5px;}
      tr:nth-child(even) {background-color: #dddddd;}
      </style> </head> <body>
      <table>";
    for(j in 1:nrow(transposed_df)){
      text <- gsub(",", " ", transposed_df[j,i])
      text <- gsub("(.{101,}?)\\s", "\\1<br>", text)
      node_title = paste0(node_title,
                          " <tr>",
                          "<td>",
                          names_rows[j],
                          "</td>",
                          "<td>",
                          text,
                          "</td>",
                          "</tr> ")
    }
    node_title = paste0(node_title, "</table> </body> </html>")
    title[[i]] <- node_title
  }




  V(g)$title[ids] <- paste0(
    "<h4>",
    sprintf('<a href="http://amigo.geneontology.org/amigo/term/%s" target="_blank">%s</a>', gs_names[ids], gs_names[ids]), "</h4><br>",
    title[ids], "<br><br>")
  return(g)

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
#' genes <- list(c("PDHB", "VARS2"), c("IARS2", "PDHA1"),
#'  c("AAAS", "ABCE1"), c("ABI1", "AAR2"), c("AATF", "AMFR"),
#'  c("BMS1", "DAP3"), c("AURKAIP1", "CHCHD1"), c("IARS2"),
#'  c("AHI1", "ALMS1"))
#'
#' g <- getBipartiteGraph(cluster, geneset_names, genes)
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

  #cluster_size <- sapply(cluster, length)
  #igraph::V(graph)$size[cluster_id] <- 5 * sqrt(cluster_size)


  igraph::V(graph)$title <- NA

  text <- list()
  for(i in cluster_id){
    mem <- paste(geneset_names[cluster[[i]]], collapse = " ")
    mem <- gsub("(.{21,}?)\\s", "\\1<br>", mem)
    text[[i]] <- mem
  }

  for(j in geneset_id){
    gs <- paste(unlist(genes[as.integer(na.omit(rownames(df_node_mapping)[df_node_mapping$Node_number == j]))]
                ),collapse = " ")
    gs <- gsub("(.{71,}?)\\s", "\\1<br>", gs)
    text[[j]] <- gs
  }

  for(i in cluster_id){
    igraph::V(graph)$title[i] <- paste(
      "<h4>", igraph::V(graph)$name[i], "</h4><br>",
      "Members:<br>", text[[i]]
    )
  }

  for(j in geneset_id){
    igraph::V(graph)$title[[j]] <- paste(
      "<h4>", igraph::V(graph)$name[j], "</h4><br>",
      "Genes:<br>", text[[j]])
  }

  return(graph)
}
