#' Check for subset inclusion
#'
#' Remove subsets from a given list of sets, i.e. remove sets which are
#' completely contained in any other larger set in the list.
#'
#' @param seeds A `list` of sets
#'
#' @return A `list` of unique sets
#' @importFrom utils combn
#' @export
#'
#' @examples
#' ## Mock example showing how the data should look like
#'
#' seeds <- list(c(1:5), c(2:5), c(6:10))
#' s <- checkInclusion(seeds)
#'
#' ## Example using the data available in the package
#' data(scores_macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#'
#' seeds <- seedFinding(scores_macrophage_topGO_example_small,
#'                      simThreshold = 0.3,
#'                      memThreshold = 0.5)
#' seeds <- checkInclusion(seeds)
checkInclusion <- function(seeds) {
  # Remove all empty seeds from the list
  seeds <- seeds[lengths(seeds) != 0]
  
  # Create a list to store the indices of sets to be removed (i.e. subsets of
  # other sets)
  remove <- c()
  
  # If there are less than 2 sets, no removal is possible,
  # return the original set list
  if (length(seeds) < 2) {
    return(seeds)
  }
  
  # Determine the number of sets
  l <- length(seeds)
  
  # Iterate over all sets to compare them for inclusion
  for (i in seq_len((l))) {
    if (i %in% remove) {
      next
    }
    # Current set
    s1 <- seeds[[i]]
    l1 <- length(s1)
    comb <-
      lapply(
        c(1:l1),
        FUN = function(x)
          combn(s1, x, simplify = F)
      )
    comb <- do.call(c, comb)
    r <- which(seeds %in% comb)
    r <- r[r != i]
    remove <- c(remove, r)
    
  }
  # Ensure that there are no duplicates in the sets to remove
  remove <- unique(remove)
  
  # Remove the identified subsets from the original list of seeds
  if (length(remove) == 0) {
    # If there are no sets to be removed, return the original list
    return(seeds)
  } else {
    # Return the set list with subsets removed
    return(seeds[-remove])
  }
}


#' Find clustering seeds
#'
#' Determine initial seeds for the clustering from the distance score matrix.
#'
#' @references
#' See https://david.ncifcrf.gov/helps/functional_classification.html#clustering
#' for details on the original implementation
#'
#' @param distances A [Matrix::Matrix()] of (distance) scores
#' @param simThreshold numerical, A threshold to determine which genesets are
#'                     considered close (i.e. have a distance <= simThreshold)
#'                     in the `distances` matrix.
#' @param memThreshold numerical, A threshold used to ensure that enough members
#'                     of a potential seed set are close/similar to each other.
#'                     Only if this condition is met, the set is considered a
#'                     seed.
#'
#' @return A `list` of seeds which can be used for clustering
#' @export
#'
#' @examples
#' ## Mock example showing how the data should look like
#'
#' m <- Matrix::Matrix(stats::runif(100, min = 0, max = 1), 10, 10)
#' seeds <- seedFinding(distances = m, simThreshold = 0.3, memThreshold = 0.5)
#'
#' ## Example using the data available in the package
#' data(scores_macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#'
#' seeds <- seedFinding(scores_macrophage_topGO_example_small,
#'                      simThreshold = 0.3,
#'                      memThreshold = 0.5)
seedFinding <- function(distances, simThreshold, memThreshold) {
  # Check if there are any distance scores, if not, return NULL
  if (is.null(distances) || length(distances) == 0) {
    return(NULL)
  }
  
  # Initialize a list to store the identified seeds
  seeds <- list()
  
  # Determine which entries of the distances matrix are reachable from each
  # other (i.e. have a distance score smaller or equal the provided
  # simThreshold)
  reach <-
    apply(distances, 1, function(x) {
      as.numeric(x <= simThreshold)
    })
  
  # Iterate over all rows in the distance score matrix
  for (i in seq_len(nrow(distances))) {
    # Check if at least 2 other entries are reachable from i
    if (sum(reach[i, ], na.rm = TRUE) >= 2) {
      # Extract members reachable from i
      members <- which(reach[i, ] == 1)
      
      # Calculate an individual threshold for i to be considered a seed
      includethreshold <-
        (length(members) ^ 2 - length(members)) * memThreshold
      
      # Subset the reach matrix and sum up entries
      reach_red <- reach[members, members]
      in_reach <- sum(reach_red)
      
      # If sum of entries in reach is above the individual threshold,
      # i is a seed
      if (in_reach >= includethreshold) {
        members <- c(members, i)
        seeds <- c(list(sort(members)), seeds)
      }
    }
  }
  
  seeds <- unique(seeds)
  # Ensure each seed contains each member only once
  seeds <- lapply(seeds, unique)
  # Ensure no seed is fully included in a larger seed
  seeds <- checkInclusion(seeds)
  
  # Return the identified seeds
  return(seeds)
}


#' Find cluster from initial seeds
#'
#' Merge the initially determined seeds to clusters.
#'
#' @references
#' See https://david.ncifcrf.gov/helps/functional_classification.html#clustering
#' for details on the original implementation
#'
#' @param seeds A `list` of seeds, e.g. determined by \code{GeDi::seedFinding()}
#'              function
#' @param threshold numerical, A threshold for merging seeds
#'
#' @return A `list` of clusters
#' @export
#'
#' @examples
#' ## Mock example showing how the data should look like
#'
#' seeds <- list(c(1:5), c(6:10))
#' cluster <- fuzzyClustering(seeds, 0.5)
#'
#' ## Example using the data available in the package
#' data(scores_macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#'
#' seeds <- seedFinding(scores_macrophage_topGO_example_small,
#'                      simThreshold = 0.3,
#'                      memThreshold = 0.5)
#' cluster <- fuzzyClustering(seeds, threshold = 0.5)
fuzzyClustering <- function(seeds, threshold) {
  # Check if there are at least two seeds to merge
  # If not, return the original seeds
  if (length(seeds) <= 1) {
    return(seeds)
  }
  
  # Create a logical vector to track whether a seed is still mergeable
  mergeable <- rep(TRUE, length(seeds))
  
  # Repeat the merging process until no more seeds are mergeable
  while (any(mergeable)) {
    # Get the index of the first mergeable seed
    index <- which(mergeable)[1]
    if (index > length(seeds)) {
      break
    }
    
    # Current mergeable seed
    s1 <- seeds[[index]]
    l <- length(seeds)
    
    # Iterate over all seeds to check for merging possibilities
    for (j in seq_len(length(seeds))) {
      s2 <- seeds[[j]]
      int <- intersect(s1, s2)
      union <- sort(union(s1, s2))
      # Check if the two seeds are mergeable
      if (length(int) >= (threshold * length(union))) {
        # If mergeable, remove the individual seeds and add a new merged seed
        remove <- list(s1, s2)
        seeds <- seeds[!(seeds %in% remove)]
        seeds <- c(list(union), seeds)
        mergeable <- mergeable[-c(index, j)]
        mergeable <- c(TRUE, mergeable)
        break
      }
    }
    # Check if there are still seeds to merge, otherwise mark the seed
    # as unmergeable
    if (l == length(seeds)) {
      mergeable[[index]] <- FALSE
    }
  }
  
  return(seeds)
}


#' Cluster genesets.
#'
#' This function performs clustering on a set of scores using either the Louvain
#' or Markov method.
#'
#' @param scores A [Matrix::Matrix()] of (distance) scores
#' @param threshold numerical, A threshold used to determine which genesets are
#'                  considered similar. Genesets are considered similar if
#'                  (distance) score <= threshold.
#'                  similar.
#' @param cluster_method character, the clustering method to use. The options
#'                       are `louvain` and `markov`. Defaults to `louvain`.
#'
#' @return A `list` of clusters
#' @export
#' @importFrom igraph cluster_louvain membership
#' @importFrom GeneTonic cluster_markov
#'
#' @examples
#' ## Mock example showing how the data should look like
#' m <- Matrix::Matrix(stats::runif(100, min = 0, max = 1), 10, 10)
#' rownames(m) <- colnames(m) <- c("a", "b", "c", "d", "e",
#'                                 "f", "g", "h", "i", "j")
#' cluster <- clustering(m, 0.3, "markov")
#'
#' ## Example using the data available in the package
#' data(scores_macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#'
#'clustering <- clustering(scores_macrophage_topGO_example_small,
#'                         threshold = 0.5)
clustering <-
  function(scores, threshold, cluster_method = "louvain") {
    # Check if the cluster_method is valid (only "louvain" or "markov" allowed)
    stopifnot(cluster_method == "louvain" ||
                cluster_method == "markov")
    
    # Obtain adjacency matrix based on the distance scores and build a graph
    adj_matrix <- getAdjacencyMatrix(scores, threshold)
    stopifnot(!is.null(adj_matrix))
    graph <- buildGraph(adj_matrix)
    
    # Run Louvain or Markov clustering based on the chosen method
    if (cluster_method == "louvain") {
      clustering <- cluster_louvain(graph)
      memberships <- membership(clustering)
    } else if (cluster_method == "markov") {
      clustering <- cluster_markov(graph)
      memberships <- clustering$membership
    }
    
    # Extract cluster memberships for each geneset
    cluster <- vector(mode = "list", length = max(memberships))
    
    # Transform the mapping of geneset -> cluster to cluster -> genesets mapping
    for (i in seq_len(length(memberships))) {
      sub_cluster <- memberships[i]
      cluster[[sub_cluster]] <- c(cluster[[sub_cluster]], i)
    }
    
    # Remove all singleton clusters (clusters with only one geneset)
    filter <- vapply(cluster, function(x)
      length(x) > 1, logical(1))
    cluster <- cluster[filter]
    
    # Return the final cluster mapping
    return(cluster)
  }


#' Calculate clusters based on kNN clustering
#'
#' This function performs k-Nearest Neighbors (kNN) clustering on a set of
#' scores.
#'
#' @param scores A [Matrix::Matrix()] of (distance) scores
#' @param k numerical, the number of neighbors
#'
#' @return A `list` of clusters
#' @export
#' @importFrom BiocNeighbors findKNN
#'
#' @examples
#' ## Mock example showing how the data should look like
#' scores <- Matrix::Matrix(stats::runif(100, min = 0, max = 1), 10, 10)
#' rownames(scores) <- colnames(scores) <- c("a", "b", "c", "d", "e",
#'                                 "f", "g", "h", "i", "j")
#' cluster <- kNN_clustering(scores, k = 3)
#'
#' ## Example using the data available in the package
#' data(scores_macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#'
#'kNN <- kNN_clustering(scores_macrophage_topGO_example_small,
#'                         k = 5)
kNN_clustering <- function(scores,
                           k) {
  # Check if there are any distance scores, if not, return NULL
  if (is.null(scores) || length(scores) == 0) {
    return(NULL)
  }
  # Find k nearest neighbors for each geneset in the data
  kNN <- findKNN(scores, k, warn.ties = FALSE)
  # Extract the list of neighbors for each geneset
  kNN <- kNN$index
  # Select the first neighbor as the cluster for each geneset
  kNN <- lapply(seq_len(nrow(kNN)), function(i)
    kNN[i, ])
  
  # Return the list of clusters based on k-Nearest Neighbors
  return(kNN)
}


#' Map each geneset to the cluster it belongs
#'
#' Map each geneset to the cluster it belongs and return the information as
#' a `data.frame`
#'
#' @param cluster A `list` of clusters
#' @param gs_names A vector of geneset names
#' @param gs_description A vector of descriptions for each geneset
#'
#' @return A `data.frame` mapping each geneset to the cluster(s) it belongs to
.getClusterDatatable <-
  function(cluster, gs_names, gs_description) {
    # Check if geneset names are given
    stopifnot(length(gs_names) > 0)
    n_gs <- length(gs_names)
    df <- vector("list", n_gs)
    
    # Check if there are no clusters
    if (length(cluster) == 0) {
      # Create a data.frame with "No associated Cluster" label for all genesets
      df <-
        data.frame(Cluster = rep("No associated Cluster", n_gs))
      rownames(df) <- gs_names
      return(df)
    }
    
    # Iterate over all clusters and genesets to build up the data.frame
    for (i in seq_len(length(cluster))) {
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
    
    # Transform the list of lists into a data.frame format
    df <- data.frame(matrix(df, nrow = n_gs, ncol = 1))
    colnames(df) <- c("Cluster")
    cluster <- df$Cluster
    
    # Set information of genesets belonging to no cluster
    cluster <- lapply(cluster, function(x) {
      if (is.null(x)) {
        "No associated Cluster"
      } else {
        x
      }
    })
    
    df$Cluster <- cluster
    df$Description <- gs_description
    df <- df[, c(2, 1)]
    rownames(df) <- gs_names
    
    # Return the final cluster datatable
    return(df)
  }
