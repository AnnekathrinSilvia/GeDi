#' Check for subset inclusion
#'
#' Check if every set is unique and there are no subset included in larger sets
#'
#' @param seeds A `list` of sets of numerical values
#'
#' @return A `list` of unique sets
#' @export
#'
#' @examples
#' seeds <- list(c(1:5), c(2:5), c(6:10))
#' s <- checkInclusion(seeds)
checkInclusion <- function(seeds) {
  # set up a list for the sets to remove (i.e. those completely
  # included in other sets)

  remove <- c()

  # if there are not at least 2 sets, nothing can be removed
  if (length(seeds) < 2) {
    return(seeds)
  }

  l <- length(seeds)

  # iterate over all sets
  for (i in 1:(l - 1)) {
    s1 <- seeds[[i]]
    # iterate over all remaining sets
    for (j in (i + 1):l) {
      s2 <- seeds[[j]]

      # check if either s1 or s2 is a subset of the other, if so add to the
      # list of sets to remove
      if (setequal(intersect(s1, s2), s1)) {
        remove <- c(remove, i)
      } else if (setequal(intersect(s2, s1), s2)) {
        remove <- c(remove, j)
      }
    }
  }

  # check that no set is duplicated in remove
  remove <- unique(remove)

  # remove the sets from the initial list
  if (length(remove) == 0) {
    return(seeds)
  } else {
    return(seeds[-remove])
  }
}

#' Find clustering seeds
#'
#' Find initial seeds for the clustering. #TODO: Verlinken auf die David publication
#'
#' @param distances A [Matrix::Matrix()] of (distance) scores
#' @param simThreshold numerical, a threshold of what is considered a
#'                     close relationship between two elements in distances
#' @param memThreshold numerical, a threshold used to identify members of
#'                     a seeds
#'
#' @return A `list` of seeds which can be used for clustering
#' @export
#'
#' @examples
#' m <- Matrix::Matrix(stats::runif(100, min = 0, max = 1), 10, 10)
#' seeds <- seedFinding(distances = m, simThreshold = 0.3, memThreshold = 0.5)
seedFinding <- function(distances, simThreshold, memThreshold) {
  # first check if there a any distance scores
  if (is.null(distances) || length(distances) == 0) {
    return(NULL)
  }

  # set up a list for the seeds
  seeds <- list()

  # check which entries of distances are reachable from each other
  # (i.e. have a distance score smaller/eqaul simThreshold)
  reach <-
    apply(distances, 1, function(x) {
      as.numeric(x <= simThreshold)
    })

  # iterate over all rows in the distance score
  for (i in 1:nrow(distances)) {
    # check if at least 2 other entries are reachable from i
    if (sum(reach[i, ], na.rm = TRUE) >= 2) {
      # extract the members which are reachable from i
      members <- which(reach[i, ] == 1)
      # calculate the individual threshold for i to be considered a seed
      includethreshold <-
        (length(members)^2 - length(members)) * memThreshold
      # subset the reach matrix and sum up the entries in reach
      reach_red <- reach[members, members]
      in_reach <- sum(reach_red)
      # if the sum of entries in reach is larger than the individual threshold
      # of i the set is considered a seed
      if (in_reach >= includethreshold) {
        members <- c(members, i)
        seeds <- c(list(sort(members)), seeds)
      }
    }
  }

  # check than no seed is completely included in a larger seed
  seeds <- checkInclusion(unique(seeds))
  # check that each seed contains each member only once
  seeds <- lapply(seeds, unique)
  return(seeds)
}


#' Find cluster from initial seeds
#'
#' Find clusters from initial seeds according to (#TODO: David publication linken)
#'
#' @param seeds A `list` of seeds, e.g. determined by seedFinding
#' @param threshold numerical, threshold for merging seeds
#'
#' @return A `list` of clusters
#' @export
#'
#' @examples
#' seeds <- list(c(1:5), c(6:10))
#' cluster <- fuzzy_clustering(seeds, 0.5)
fuzzy_clustering <- function(seeds, threshold) {
  # Check if there seeds to merge
  if (length(seeds) <= 1) {
    return(seeds)
  }

  # set a logical vector to check if a seed is still mergeable
  mergeable <- rep(TRUE, length(seeds))

  # repeat until there are no longer mergeable seeds
  while (any(mergeable)) {
    # get the first mergeable seed
    index <- which(mergeable)[1]
    if (index > length(seeds)) {
      break
    }
    s1 <- seeds[[index]]
    l <- length(seeds)
    for (j in 1:length(seeds)) {
      s2 <- seeds[[j]]
      int <- intersect(s1, s2)
      union <- sort(union(s1, s2))
      # check if the two seeds are mergeable according to a specific condition
      if (length(int) >= (threshold * length(union))) {
        # if mergeable, remove the two individual seeds from the list of seeds
        # and add a new merged seed
        remove <- list(s1, s2)
        seeds <- seeds[!(seeds %in% remove)]
        seeds <- c(list(union), seeds)
        mergeable <- mergeable[-c(index, j)]
        mergeable <- c(TRUE, mergeable)
        break
      }
    }
    # check if there are still seeds to merge
    if (l == length(seeds)) {
      mergeable[[index]] <- FALSE
    }
  }

  return(seeds)
}

#' Cluster genesets using Louvain or Markov clustering.
#'
#' Cluster the geneset using either Louvain or Markov clustering.
#'
#' @param scores A [Matrix::Matrix()] of (distance) scores
#' @param threshold numerical, a threshold indicating similar genesets. Genesets
#'                  with a (distance) score <= threshold will be considered
#'                  similar.
#' @param cluster_method character, the clustering method to use. The options
#'                       are `louvain` and `markov`. Default to `louvain`.
#'
#' @return A `list` of clusters
#' @export
#' @importFrom igraph cluster_louvain membership
#' @importFrom GeneTonic cluster_markov
#'
#' @examples
#' m <- Matrix::Matrix(stats::runif(100, min = 0, max = 1), 10, 10)
#' rownames(m) <- colnames(m) <- c("a", "b", "c", "d", "e",
#'                                 "f", "g", "h", "i", "j")
#' cluster <- clustering(m, 0.3, "markov")
clustering <- function(scores,
                       threshold,
                       cluster_method = "louvain"){
  stopifnot(cluster_method == "louvain" || cluster_method == "markov")
  # get adjacency matrix of the data and build a graph
  adj_matrix <- getAdjacencyMatrix(scores,
                                   threshold)
  graph <- buildGraph(adj_matrix)

  # run louvain or markov clustering
  if(cluster_method == "louvain"){
    clustering <- cluster_louvain(graph)
    memberships <- membership(clustering)
  }else if(cluster_method == "markov"){
    clustering <- cluster_markov(graph)
    memberships <- clustering$membership
  }

  # extract the cluster memeberships of each genesets
  cluster <- vector(mode = "list", length = max(memberships))

  # tranform the mapping of geneset -> cluster to a cluster-> genesets mapping
  for(i in 1:length(memberships)){
    sub_cluster <- memberships[i]
    cluster[[sub_cluster]] <- c(cluster[[sub_cluster]], i)
  }

  # remove all singletons
  filter <- sapply(cluster, function(x) length(x) > 1)
  cluster <- cluster[filter]

  return(cluster)
}

#' Calculate clusters based on kNN clustering
#'
#' Calculate a clustering of the data using the kNN approach
#'
#' @param scores A [Matrix::Matrix()] of (distance) scores
#' @param k numerical, the number of neighbors that should be searched for
#'
#' @return A `list` of clusters
#' @export
#' @import BiocNeighbors
#'
#' @examples
#' scores <- Matrix::Matrix(stats::runif(100, min = 0, max = 1), 10, 10)
#' rownames(scores) <- colnames(scores) <- c("a", "b", "c", "d", "e",
#'                                 "f", "g", "h", "i", "j")
#' cluster <- kNN_clustering(scores, k = 3)
kNN_clustering <- function(scores,
                           k){
  # find k nearest neighbors for each geneset in the data
  kNN <- findKNN(scores, k)
  # extract for each geneset the list of neighbors
  kNN <- kNN$index
  # select the first neighbor as cluster for each geneset
  kNN <- lapply(seq_len(nrow(kNN)), function(i) kNN[i,])
  return(kNN)
}


#' Map each geneset to the cluster it belongs
#'
#' Map each geneset to the cluster it belongs and return the information as
#' a `data.frame`
#'
#' @param cluster A `list` of clusters
#' @param gs_names A vector of geneset names
#'
#' @return A `data.frame` mapping each geneset to the cluster(s) it belongs to
.getClusterDatatable <- function(cluster, gs_names) {
  # check if geneset names are given
  stopifnot(length(gs_names) > 0)
  n_gs <- length(gs_names)
  df <- vector("list", n_gs)

  # check if there are any clusters
  if (length(cluster) == 0) {
    df <- data.frame(Cluster = rep("No associated Cluster", n_gs))
    rownames(df) <- gs_names
    return(df)
  }

  # iterate over all clusters and build up the data.frame
  for (i in 1:length(cluster)) {
    for (j in cluster[[i]]) {
      entry <- df[[j]]
      if (is.null(entry)) {
        entry <- i
      } else {
        entry <- c(entry, i)
      }

      df[[j]] <- entry
    }
  }

  # transfrom into desired data.frame format
  df <- data.frame(matrix(df, nrow = n_gs, ncol = 1))
  colnames(df) <- c("Cluster")
  cluster <- df$Cluster
  # set information of genesets belonging to no cluster
  cluster <-
    lapply(cluster, function(x) {
      if (is.null(x)) {
        "No associated Cluster"
      } else {
        x
      }
    })
  df$Cluster <- cluster
  rownames(df) <- gs_names
  return(df)
}
