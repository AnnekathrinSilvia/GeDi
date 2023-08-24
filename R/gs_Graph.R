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
#' m <- Matrix::Matrix(stats::runif(1000, 0, 1), 100, 100)
#' geneset_names <- as.character(stats::runif(100, min = 0, max = 1))
#' rownames(m) <- colnames(m) <- geneset_name
#' threshold <- 0.3
#' adj <- getAdjacencyMatrix(m, threshold)
getAdjacencyMatrix <- function(distanceMatrix,
                               cutOff) {
  # Ensure that the distance score matrix is valid
  if (is.null(distanceMatrix)) {
    return(NULL)
  }

  # Determine the number of nodes, which is equal to the number of rows
  l <- nrow(distanceMatrix)
  # Initialize an adjacency matrix with zeros
  adjMat <- Matrix(0, l, l)

  # Iterate over each node to identify edges based on cutoff
  for (i in 1:l) {
    # Get indices of nodes within cutoff
    edge <- which(distanceMatrix[i,] <= cutOff)
    if (length(edge) > 0) {
      # Set adjacency matrix entry to 1 for connected nodes
      adjMat[i, edge] <- 1
    }
  }

  # Set row and column names of the adjacency matrix
  rownames(adjMat) <- rownames(distanceMatrix)
  colnames(adjMat) <- colnames(distanceMatrix)

  # Return the adjacency matrix
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
#' adj <- Matrix::Matrix(0, 100, 100)
#' adj[c(80:100), c(80:100)] <- 1
#' geneset_names <- as.character(stats::runif(100, min = 0, max = 1))
#' rownames(adj) <- colnames(adj) <- geneset_names
#' graph <- buildGraph(adj)
buildGraph <- function(adjMatrix) {
  # Build an undirected graph from the adjacency matrix
  g <- graph_from_adjacency_matrix(
    adjMatrix,
    mode = "undirected",
    add.colnames = NULL,
    add.rownames = NA,
    diag = F
  )

  # Get geneset names from row names of adjacency matrix
  gs_names <- rownames(adjMatrix)
  # Get indices of nodes that match geneset names
  ids <- which(names(V(g)) %in% gs_names)

  # Customize node titles based on the type of database (GO, Reactome, or other)
  if (all(sapply(gs_names, function(x)
    substr(x, 1, 2) == "GO"))) {
    V(g)$title[ids] <- paste0(
      "<h4>",
      sprintf(
        '<a href="http://amigo.geneontology.org/amigo/term/%s" target="_blank">%s</a>',
        gs_names[ids],
        gs_names[ids]
      ),
      "</h4><br>",
      V(g)$name[ids],
      "<br><br>"
    )
  } else if (all(sapply(gs_names, function(x)
    substr(x, 1, 2) == "R-"))) {
    V(g)$title[ids] <- paste0(
      "<h4>",
      sprintf(
        '<a href="http://reactome.org/content/detail/%s" target="_blank">%s</a>',
        gs_names[ids],
        gs_names[ids]
      ),
      "</h4><br>",
      V(g)$name[ids],
      "<br><br>"
    )
  } else{
    V(g)$title[ids] <- paste0(
      "<h4>",
      sprintf(
        '<a href="https://www.genome.jp/dbget-bin/www_bget?pathway:%s" target="_blank">%s</a>',
        gs_names[ids],
        gs_names[ids]
      ),
      "</h4><br>",
      V(g)$name[ids],
      "<br><br>"
    )
  }

  # Return the customized graph
  return(g)
}

#' Construct an adjacency matrix
#'
#' Construct an adjacency matrix from a `list` of cluster.
#'
#' @param cluster A `list` of clusters, where each cluster member is indicated
#'                by a numeric value
#' @param gs_names A vector of geneset names
#'
#' @return A [Matrix::Matrix()] of adjacency status
#' @importFrom Matrix Matrix
#' @export
#'
#' @examples
#' cluster <- list(c(1:5), c(6:9))
#' gs_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
#' adj <- getClusterAdjacencyMatrix(cluster, gs_names)
getClusterAdjacencyMatrix <- function(cluster,
                                      gs_names) {
  # Number of genesets
  l <- length(gs_names)
  # Ensure that there is at least one geneset
  stopifnot(l > 0)

  # Initialize an adjacency matrix with zeros
  adj <- Matrix::Matrix(0, l, l)

  # Check if the cluster is empty
  if (length(cluster) == 0) {
    rownames(adj) <- colnames(adj) <- gs_names
    diag(adj) <- 0
    return(adj)
  }

  # Ensure that there are not more cluster than genesets
  stopifnot(l >= max(unlist(cluster)))

  # Fill the adjacency matrix based on the provided cluster
  for (i in 1:length(cluster)) {
    # Get subcluster indices
    subcluster <- cluster[[i]]
    # Initialize edges between all nodes in a cluster
    adj[subcluster, subcluster] <- 1
  }

  # Set row and column names
  rownames(adj) <- colnames(adj) <- gs_names
  # Remove self-loops
  diag(adj) <- 0

  # Return the adjacency matrix
  return(adj)
}

#' Build a cluster graph
#'
#' Build a [igraph] from cluster information, connecting nodes which belong to
#' the same cluster.
#'
#' @param cluster list, a `list` of clusters, where each cluster member is
#'                indicated by a numeric value.
#' @param geneset_df `data.frame`, a `data.frame` of genesets with at least two
#'                   columns, one called `Genesets` containing geneset
#'                   identifiers and one called `Genes` containing a list of
#'                   genes belonging to the individual genesets.
#' @param gs_names vector, a vector of geneset identifiers, e.g. the `Genesets`
#'                 column of `geneset_df`.
#' @param color_by character, a column name of `geneset_df` which is used
#'                 to color the nodes of the resulting graph. The column should
#'                 ideally contain a numeric measurement. Defaults to NULL and
#'                 nodes will remain uncolored.
#'
#' @return An `igraph` object to be further manipulated or processed/plotted
#'         (e.g. via [igraph::plot.igraph()] or
#'         [visNetwork::visIgraph()][visNetwork::visNetwork-igraph])
#' @import igraph
#' @importFrom GeneTonic map2color
#' @importFrom grDevices colorRampPalette
#' @export
#'
#' @examples
#' cluster <- list(c(1:5), c(6:9, 1))
#' genes <- list(
#'   c("PDHB", "VARS2"), c("IARS2", "PDHA1"),
#'   c("AAAS", "ABCE1"), c("ABI1", "AAR2"), c("AATF", "AMFR"),
#'   c("BMS1", "DAP3"), c("AURKAIP1", "CHCHD1"), c("IARS2"),
#'   c("AHI1", "ALMS1")
#' )
#' gs_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
#' geneset_df <- data.frame(
#'   Genesets = gs_names,
#'   value = rep(1, 9)
#' )
#' geneset_df$Genes <- genes
#' graph <- buildClusterGraph(
#'   cluster = cluster,
#'   geneset_df = geneset_df,
#'   gs_names = gs_names,
#'   color_by = "value"
#' )
buildClusterGraph <- function(cluster,
                              geneset_df,
                              gs_names,
                              color_by = NULL) {
  # Get adjacency matrix representing genesets belonging to the same cluster
  adj <- getClusterAdjacencyMatrix(cluster,
                                   gs_names)
  # Build a graph from the adjacency matrix
  g <- buildGraph(adj)

  # Get node ids corresponding to geneset names
  ids <- which(names(V(g)) %in% gs_names)

  # Add cluster information to nodes in the graph
  V(g)$cluster <- ""

  n_cluster <- length(cluster)
  if (n_cluster > 0) {
    for (i in 1:n_cluster) {
      clus <- cluster[[i]]
      for (y in 1:length(clus)) {
        gs_name <- gs_names[clus[[y]]]
        id <- which(names(V(g)) %in% gs_name)
        mem <- V(g)$cluster[id]
        cluster_name <- paste("Cluster ", i, sep = "")
        if (mem != "") {
          mem <- paste(mem, cluster_name, sep = ",")
        } else {
          mem <- cluster_name
        }
        V(g)$cluster[id] <- mem
      }
    }
  }

  # Construct HTML-based title for each node using input information from geneset_df
  transposed_df <- as.data.frame(t(geneset_df))
  title <- list()
  names_rows <- rownames(transposed_df)

  for (i in 1:ncol(transposed_df)) {
    node_title <- "<!DOCTYPE html> <html> <head> <style>
      table {font-family: arial, sans-serif; font-size: 10px; border-collapse: collapse;width: 100%;} td,
      th { border: 1px solid #dddddd; text-align: center; padding: 5px;}
      tr:nth-child(even) {background-color: #dddddd;}
      </style> </head> <body>
      <table>"
    for (j in 1:nrow(transposed_df)) {
      text <- gsub(",", " ", transposed_df[j, i])
      text <- gsub("(.{101,}?)\\s", "\\1<br>", text)
      node_title <- paste0(node_title,
                           " <tr>",
                           "<td>",
                           names_rows[j],
                           "</td>",
                           "<td>",
                           text,
                           "</td>",
                           "</tr> ")
    }
    node_title <- paste0(node_title, "</table> </body> </html>")
    title[[i]] <- node_title
  }

  # Set node titles using geneset names, corresponding links, and constructed HTML titles
  V(g)$title[ids] <- paste0(
    "<h4>",
    sprintf(
      '<a href="http://amigo.geneontology.org/amigo/term/%s" target="_blank">%s</a>',
      gs_names[ids],
      gs_names[ids]
    ),
    "</h4><br>",
    title[ids],
    "<br><br>"
  )

  # Remove nodes without any connections (degree equals 0)
  no_cluster <- V(g)[degree(g) == 0]
  g <- delete_vertices(g, no_cluster)

  # Update ids to include only nodes present in the graph
  ids <- which(names(V(g)) %in% gs_names)

  if (!is.null(color_by)) {
    # Check if the specified color_by column exists in geneset_df
    if (!color_by %in% colnames(geneset_df)) {
      stop(
        "Your data does not contain the ",
        color_by,
        " column.\n",
        "Please select another column to use for the color."
      )
    }
    col_var <- geneset_df[ids, color_by]
    # the palette changes if it is z_score VS pvalue
    if (all(col_var <= 1) & all(col_var > 0)) {
      # likely p-values...
      col_var <- -log10(col_var)
      mypal <- (scales::alpha(colorRampPalette(
        RColorBrewer::brewer.pal(name = "YlOrRd", 7)
      )(10), 1))
      mypal_hover <- (scales::alpha(colorRampPalette(
        RColorBrewer::brewer.pal(name = "YlOrRd", 7)
      )(10), 0.5))
      mypal_select <- (scales::alpha(colorRampPalette(
        RColorBrewer::brewer.pal(name = "YlOrRd", 7)
      )(10), 1))

      V(g)$color.background <- map2color(col_var,
                                         mypal,
                                         symmetric = FALSE,
                                         limits = range(na.omit(col_var)))
      V(g)$color.highlight <- map2color(col_var,
                                        mypal_select,
                                        symmetric = FALSE,
                                        limits = range(na.omit(col_var)))
      V(g)$color.hover <- map2color(col_var,
                                    mypal_hover,
                                    symmetric = FALSE,
                                    limits = range(na.omit(col_var)))

      V(g)$color.background[is.na(V(g)$color.background)] <-
        "lightgrey"
        V(g)$color.highlight[is.na(V(g)$color.highlight)] <-
          "lightgrey"
          V(g)$color.hover[is.na(V(g)$color.hover)] <- "lightgrey"
    } else {
      # e.g. using z_score or aggregated value
      if (prod(range(na.omit(col_var))) >= 0) {
        # gradient palette
        mypal <- (scales::alpha(colorRampPalette(
          RColorBrewer::brewer.pal(name = "Reds", 5)
        )(5), 1))
        mypal_hover <- (scales::alpha(colorRampPalette(
          RColorBrewer::brewer.pal(name = "Reds", 5)
        )(5), 0.5))
        mypal_select <- (scales::alpha(colorRampPalette(
          RColorBrewer::brewer.pal(name = "Reds", 5)
        )(5), 1))

        V(g)$color.background <- map2color(col_var,
                                           mypal,
                                           symmetric = FALSE,
                                           limits = range(na.omit(col_var)))
        V(g)$color.highlight <- map2color(col_var,
                                          mypal_select,
                                          symmetric = FALSE,
                                          limits = range(na.omit(col_var)))
        V(g)$color.hover <- map2color(col_var,
                                      mypal_hover,
                                      symmetric = FALSE,
                                      limits = range(na.omit(col_var)))
        V(g)$color.background[is.na(V(g)$color.background)] <-
          "lightgrey"
          V(g)$color.highlight[is.na(V(g)$color.highlight)] <-
            "lightgrey"
            V(g)$color.hover[is.na(V(g)$color.hover)] <- "lightgrey"
      }
      #  else {
      #   # divergent palette to be used
      #   mypal <- rev(scales::alpha(
      #     colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.8
      #   ))
      #   mypal_hover <- rev(scales::alpha(
      #     colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.5
      #   ))
      #   mypal_select <- rev(scales::alpha(
      #     colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 1
      #   ))
      #
      #   V(g)$color.background <- map2color(col_var, mypal,
      #     symmetric = TRUE,
      #     limits = range(na.omit(col_var))
      #   )
      #   V(g)$color.highlight <- map2color(col_var, mypal_select,
      #     symmetric = TRUE,
      #     limits = range(na.omit(col_var))
      #   )
      #   V(g)$color.hover <- map2color(col_var, mypal_hover,
      #     symmetric = TRUE,
      #     limits = range(na.omit(col_var))
      #   )
      #
      #   V(g)$color.background[is.na(V(g)$color.background)] <- "lightgrey"
      #   V(g)$color.highlight[is.na(V(g)$color.highlight)] <- "lightgrey"
      #   V(g)$color.hover[is.na(V(g)$color.hover)] <- "lightgrey"
      # }
    }
  }

  # Return the constructed graph
  return(g)
}

#' Construct a bipartite graph
#'
#' Construct a bipartite graph from cluster information, mapping the cluster
#' to its members
#'
#' @param cluster `list`, a `list` of clusters, cluster members are indicated by
#'                numeric values.
#' @param gs_names vector, a vector of (geneset) identifiers/names to map the
#'                 numeric member value in `cluster` to.
#' @param genes `list`, a `list` of vectors of genenames which belong to the
#'               genesets in `gs_names`.
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
#' gs_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
#' genes <- list(
#'   c("PDHB", "VARS2"), c("IARS2", "PDHA1"),
#'   c("AAAS", "ABCE1"), c("ABI1", "AAR2"), c("AATF", "AMFR"),
#'   c("BMS1", "DAP3"), c("AURKAIP1", "CHCHD1"), c("IARS2"),
#'   c("AHI1", "ALMS1")
#' )
#'
#' g <- getBipartiteGraph(cluster, gs_names, genes)
getBipartiteGraph <- function(cluster,
                              gs_names,
                              genes) {
  # Check if parameters fulfill requirements
  stopifnot(length(cluster) > 0)
  stopifnot(length(gs_names) > 0)
  stopifnot(length(genes) > 0)

  edgelist <- c()
  type <- c()
  n_cluster <- length(cluster)
  node_number <- n_cluster + 1
  df_node_mapping <-
    data.frame(matrix(NA, nrow = length(gs_names), ncol = 1))
  colnames(df_node_mapping) <- "Node_number"

  # Set up node labels for the clusters
  node_labels <- c()
  for (i in 1:n_cluster) {
    node_labels <- c(node_labels, paste0("Cluster ", i))
  }

  # Map numerical cluster member value to the respective geneset identifier
  # Add edges from cluster to members of the cluster
  for (i in 1:n_cluster) {
    subcluster <- cluster[[i]]
    for (j in subcluster) {
      if (!is.na(df_node_mapping[j,])) {
        n <- df_node_mapping[j,]
        edgelist <- c(edgelist, i, n)
      } else {
        edgelist <- c(edgelist, i, node_number)
        df_node_mapping[j,] <- node_number
        node_number <- node_number + 1
        node_labels <- c(node_labels, gs_names[[j]])
      }
    }
  }

  # Set up the bipartite graph
  type <- c(rep(0, n_cluster))
  type <- c(type, rep(1, node_number - n_cluster - 1))

  graph <-
    igraph::make_bipartite_graph(type, edgelist, directed = TRUE)
  graph <- set_vertex_attr(graph, "name", value = node_labels)

  cluster_id <- which(names(V(graph)) %in% node_labels[1:n_cluster])
  geneset_id <-
    which(!(names(V(graph)) %in% node_labels[1:n_cluster]))

  # Set node type and shape attributes
  igraph::V(graph)$nodeType <- NA
  igraph::V(graph)$nodeType[cluster_id] <- "Cluster"
  igraph::V(graph)$nodeType[geneset_id] <- "Geneset"

  igraph::V(graph)$shape <-
    c("box", "ellipse")[factor(V(graph)$nodeType, levels = c("Cluster", "Geneset"))]

  # Set color attributes for nodes and edges
  igraph::V(graph)$color <- NA
  igraph::V(graph)$color[cluster_id] <- "gold"
  igraph::V(graph)$color[geneset_id] <- "#0092AC"
  igraph::E(graph)$color <- "black"


  # Set title information for each node
  igraph::V(graph)$title <- NA
  text <- list()
  for (i in cluster_id) {
    mem <- paste(gs_names[cluster[[i]]], collapse = " ")
    mem <- gsub("(.{21,}?)\\s", "\\1<br>", mem)
    text[[i]] <- mem
  }

  for (j in geneset_id) {
    gs <-
    paste(unlist(genes[as.integer(na.omit(rownames(df_node_mapping)[df_node_mapping$Node_number == j]))]), collapse = " ")
    gs <- gsub("(.{71,}?)\\s", "\\1<br>", gs)
    text[[j]] <- gs
  }

  for (i in cluster_id) {
    igraph::V(graph)$title[i] <- paste("<h4>",
                                        igraph::V(graph)$name[i],
                                        "</h4><br>",
                                        "Members:<br>",
                                        text[[i]])
  }

  for (j in geneset_id) {
    igraph::V(graph)$title[[j]] <- paste("<h4>",
                                          igraph::V(graph)$name[j],
                                          "</h4><br>",
                                          "Genes:<br>",
                                             text[[j]])
  }

  # Return the constructed bipartite graph
  return(graph)
}


#' Generate a `data.frame` of graph metrics
#'
#' Generate a `data.frame` of the graph metrics degree, betweenness,
#' harmonic centrality and clustering coeficient for each node
#' in a given graph.
#'
#' @param g A [igraph] graph object
#' @param genesets A `data.frame` of genesets with a column `Genesets` containing
#'                 geneset identifiers and a column `Genes` containing the
#'                 genes belonging to each geneset
#'
#' @return A `data.frame` of `geneset` extended by columns for the degree,
#'         betweenness, harmonic centrality and clustering coefficient for each
#'         geneset.
#'
.graphMetricsGenesetsDT <- function(g,
                                    genesets) {
  # Get the names of nodes (genesets) in the graph
  nodes <- igraph::V(g)$name

  # Filter and prepare geneset data
  genesets <- genesets[,!names(genesets) %in% c("Genes")]
  genesets <- genesets[genesets$Genesets %in% nodes,]

  # Compute different graph metrics using igraph functions
  clustering_coef <- igraph::transitivity(g,
                                          type = "global")
  centrality <- igraph::harmonic_centrality(g,
                                            mode = "all")
  betweenness <- igraph::betweenness(g,
                                     directed = FALSE)
  degree <- igraph::degree(g,
                           mode = "all")

  # Create a data frame to store computed metrics along with geneset information
  df <- data.frame(
    nodes,
    degree,
    round(betweenness, 2),
    centrality,
    round(clustering_coef, 2),
    genesets
  )

  # Rename columns and order the data frame by the Degree column in descending order
  rownames(df) <- NULL
  colnames(df) <- c(
    "Geneset",
    "Degree",
    "Betweenness",
    "Harmonic Centrality",
    "Clustering Coefficient",
    names(genesets)
  )
  df <- df[order(df$Degree, decreasing = TRUE),]

  # Return the computed metrics data frame
  return(df)
}
