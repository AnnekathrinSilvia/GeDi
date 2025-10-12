#' Check for subset inclusion
#'
#' Remove subsets from a given list of sets, i.e. remove sets which are
#' completely contained in any other larger set in the list.
#'
#' @param seeds A `list` of sets
#'
#' @return A `list` of unique sets
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
  for (i in seq_len((l-1))) {
    if (i %in% remove) {
      next
    }
    # Current set
    s1 <- seeds[[i]]
    l1 <- length(s1)
  # Iterate over the remaining sets to compare them with the current set
  for (j in (i + 1):l) {
    s2 <- seeds[[j]]
    # Check if s1 is a subset of s2 or vice versa, if so, mark for removal
    if (setequal(intersect(s1, s2), s1) || length(s1) == 0) {
      remove <- c(remove, i)
    } else if (setequal(intersect(s2, s1), s2) || length(s2) == 0) {
      remove <- c(remove, j)
    }
  }
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
seedFinding <- function(distances,
                        simThreshold,
                        memThreshold) {
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
fuzzyClustering <- function(seeds,
                            threshold) {
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


#' Cluster genesets using Markov clustering.
#'
#' This function is a wrapper function for the Markov clustering.
#' The actual computation of the clustering is done in the `GeDi::clustering()`
#' function. This function is mainly a wrapper function for stand-alone use of
#' GeDi functionalities to enhance user experience and allow for a clearer
#' distinction of the individual clustering algorithms.
#'
#' @param scores A [Matrix::Matrix()] of (distance) scores
#' @param threshold numerical, A threshold used to determine which genesets are
#'                  considered similar. Genesets are considered similar if
#'                  (distance) score <= threshold.
#'                  similar.
#'
#' @return A `list` of clusters
#' @export
#'
#' @examples
#' ## Mock example showing how the data should look like
#' m <- Matrix::Matrix(stats::runif(100, min = 0, max = 1), 10, 10)
#' rownames(m) <- colnames(m) <- c("a", "b", "c", "d", "e",
#'                                 "f", "g", "h", "i", "j")
#' markovCluster <- markovClustering(m, 0.3)
#'
#' ## Example using the data available in the package
#' data(scores_macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#'
#' markovCluster <- markovClustering(scores_macrophage_topGO_example_small,
#'                         threshold = 0.5)
markovClustering <- function(scores,
                             threshold){

  cluster <- clustering(scores,
                        threshold,
                        cluster_method = "markov")

  return(cluster)
}



#' Cluster genesets using Louvain clustering.
#'
#' This function is a wrapper function for the Louvain clustering.
#' The actual computation of the clustering is done in the `GeDi::clustering()`
#' function. This function is mainly a wrapper function for stand-alone use of
#' GeDi functionalities to enhance user experience and allow for a clearer
#' distinction of the individual clustering algorithms.
#'
#' @param scores A [Matrix::Matrix()] of (distance) scores
#' @param threshold numerical, A threshold used to determine which genesets are
#'                  considered similar. Genesets are considered similar if
#'                  (distance) score <= threshold.
#'                  similar.
#'
#' @return A `list` of clusters
#' @export
#'
#' @examples
#' ## Mock example showing how the data should look like
#' m <- Matrix::Matrix(stats::runif(100, min = 0, max = 1), 10, 10)
#' rownames(m) <- colnames(m) <- c("a", "b", "c", "d", "e",
#'                                 "f", "g", "h", "i", "j")
#' louvainCluster <- louvainClustering(m, 0.3)
#'
#' ## Example using the data available in the package
#' data(scores_macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#'
#' louvainCluster <- louvainClustering(scores_macrophage_topGO_example_small,
#'                         threshold = 0.5)
louvainClustering <- function(scores,
                              threshold){
  cluster <- clustering(scores,
                        threshold,
                        cluster_method = "louvain")

  return(cluster)
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
clustering <- function(scores,
                       threshold,
                       cluster_method = "louvain") {
  # Check if the cluster_method is valid (only "louvain" or "markov" allowed)
  stopifnot(cluster_method == "louvain" ||
              cluster_method == "markov")

  # Obtain adjacency matrix based on the distance scores and build a graph
  adj_matrix <- getAdjacencyMatrix(scores, threshold, weighted = TRUE)
  stopifnot(!is.null(adj_matrix))
  graph <- buildGraph(as.matrix(adj_matrix), weighted = TRUE)

  # Run Louvain or Markov clustering based on the chosen method
  if (cluster_method == "louvain") {
    clustering <- cluster_louvain(graph)
    memberships <- membership(clustering)
    } else if (cluster_method == "markov") {
      clustering <- .cluster_markov(graph)
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
  kNN <- findKNN(as.matrix(scores), k, warn.ties = FALSE)
  # Extract the list of neighbors for each geneset
  kNN <- kNN$index
  # Select the first neighbor as the cluster for each geneset
  kNN <- lapply(seq_len(nrow(kNN)), function(i)
    kNN[i, ])

  # Return the list of clusters based on k-Nearest Neighbors
  return(kNN)
}

#' Calculate clusters based on kMeans clustering
#'
#' This function performs kMeans clustering on a set of
#' scores.
#'
#' @param scores A [Matrix::Matrix()] of (distance) scores
#' @param k numerical, the number of centers to start with. This number will
#'               correlate with the resulting number of clusters.
#' @param iter numerical, number of iterations for refinement. Defaults to 500.
#' @param nstart numerical, how often the start points should be switched.
#'               Ensures a robust clustering, as clustering is influenced by the
#'               start points. Defaults to 50.
#'
#' @return A `list` of clusters
#' @importFrom stats kmeans
#' @export
#'
#' @examples
#' ## Mock example showing how the data should look like
#' scores <- Matrix::Matrix(stats::runif(100, min = 0, max = 1), 10, 10)
#' rownames(scores) <- colnames(scores) <- c("a", "b", "c", "d", "e",
#'                                 "f", "g", "h", "i", "j")
#' cluster <- kMeansClustering(scores, k = 3)
#'
#' ## Example using the data available in the package
#' data(scores_macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#'
#'cluster <- kMeansClustering(scores_macrophage_topGO_example_small,
#'                             k = 5)
kMeansClustering <- function(scores,
                              k,
                              iter = 500,
                              nstart = 50){
  # Check if there are any distance scores, if not, return NULL
  if (is.null(scores) || length(scores) == 0) {
    return(NULL)
  }
  # k has to be positive, as this will be the number of
  # resulting clusters
  stopifnot(k > 0)

  # Find k means results data
  kMeans <- kmeans(scores, k, iter, nstart)
  cluster <- c()
  cluster <- lapply(seq_len(max(unique(kMeans$cluster))),
                    function(x) c(cluster,
                                  which(kMeans$cluster == x)))
  # Return the list of clusters based on k-Nearest Neighbors
  return(cluster)
}

#' Calculate clusters based on PAM clustering
#'
#' This function performs Partioning aroung Medoids clustering on a set of
#' scores.
#'
#' @param scores A [Matrix::Matrix()] of (distance) scores
#' @param k numerical, the number of centers to start with. This number will
#'               correlate with the resulting number of clusters.
#'
#' @return A `list` of clusters
#' @importFrom cluster pam
#' @export
#'
#' @examples
#' ## Mock example showing how the data should look like
#' scores <- Matrix::Matrix(stats::runif(100, min = 0, max = 1), 10, 10)
#' rownames(scores) <- colnames(scores) <- c("a", "b", "c", "d", "e",
#'                                 "f", "g", "h", "i", "j")
#' cluster <- pamClustering(scores, k = 3)
#'
#' ## Example using the data available in the package
#' data(scores_macrophage_topGO_example_small,
#'      package = "GeDi",
#'      envir = environment())
#'
#'cluster <- pamClustering(scores_macrophage_topGO_example_small,
#'                                k = 5)
pamClustering <- function(scores,
                          k){

  # Check if there are any distance scores, if not, return NULL
  if (is.null(scores) || length(scores) == 0) {
    return(NULL)
  }
  # k has to be positive, as this will be the number of
  # resulting clusters
  stopifnot(k > 0)

  # Find pam results data
  pam <- cluster::pam(scores, k, diss = TRUE,
                      pamonce = 5)
  cluster <- c()
  cluster <- lapply(seq_len(max(unique(pam$clustering))),
                    function(x) c(cluster,
                                  which(pam$clustering == x)))
  # Return the list of clusters based on k-Nearest Neighbors
  return(cluster)
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


#' Markov Clustering (MCL) for community detection
#'
#' This function implements the Markov Clustering (MCL) algorithm for finding community
#' structure, in an analogous way to other existing algorithms in `igraph`.
#'
#' This implementation has been driven by the nice explanations provided in
#' * https://sites.cs.ucsb.edu/~xyan/classes/CS595D-2009winter/MCL_Presentation2.pdf
#' * https://medium.com/analytics-vidhya/demystifying-markov-clustering-aeb6cdabbfc7
#' * https://github.com/GuyAllard/markov_clustering (python implementation)
#'
#' More info on the MCL: https://micans.org/mcl/index.html, and
#' https://micans.org/mcl/sec_description1.html
#'
#' @references van Dongen, S.M., Graph clustering by flow simulation (2000) PhD thesis,
#' Utrecht University Repository - https://dspace.library.uu.nl/handle/1874/848
#' @references  Enright AJ, van Dongen SM, Ouzounis CA, An efficient algorithm for
#' large-scale detection of protein families (2002) Nucleic Acids Research, Volume 30,
#' Issue 7, 1 April 2002, Pages 1575–1584, https://doi.org/10.1093/nar/30.7.1575
#'
#' @param g The input graph object
#' @param add_self_loops Logical, whether to add self-loops to the matrix by
#' setting the diagonal to `loop_value`
#' @param loop_value Numeric, the value to use for self-loops
#' @param mcl_expansion Numeric, cluster expansion factor for the Markov clustering
#' iteration - defaults to 2
#' @param mcl_inflation Numeric, cluster inflation factor for the Markov clustering
#' iteration - defaults to 2
#' @param allow_singletons Logical; if `TRUE`, single isolated vertices are allowed
#' to form their own cluster. If set to `FALSE`, all clusters of size = 1 are
#' grouped in one cluster (to be interpreted as background noise).
#' @param max_iter Numeric value for the maximum number of iterations for the
#' Markov clustering
#' @param return_node_names Logical, if the graph is named and set to `TRUE`, returns
#' the node names.
#' @param return_esm Logical, controlling whether the equilibrium state matrix should be returned
#'
#' @return This function returns a `communities` object, containing the numbers of
#' the assigned membership (in the slot `membership`). Please see the
#' [igraph::communities()] manual page for additional details
#'
#' @importFrom expm "%^%"
#' @importFrom methods is
#'
#' @examples
#' library("igraph")
#' g <- make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
#' g <- add_edges(g, c(1, 6, 1, 11, 6, 11))
#' .cluster_markov(g)
#' V(g)$color <- .cluster_markov(g)$membership
#' plot(g)
.cluster_markov <- function(g,
                            add_self_loops = TRUE,
                            loop_value = 1,
                            mcl_expansion = 2,
                            mcl_inflation = 2,
                            allow_singletons = TRUE,
                            max_iter = 100,
                            return_node_names = TRUE,
                            return_esm = FALSE) {
  # g must be a graph
  if (!is(g, "igraph")) {
    stop("You need to provide an igraph object as input")
  }

  stopifnot(is.logical(add_self_loops))
  stopifnot(loop_value >= 0)
  stopifnot(mcl_expansion > 1)
  stopifnot(mcl_inflation > 1)
  stopifnot(loop_value >= 0)
  stopifnot(max_iter > 0)

  # graph to adjacency matrix
  adj_mat <- igraph::as_adjacency_matrix(g)
  adj_mat <- as.matrix(adj_mat) # to enforce non-sparse matrix

  converged <- FALSE

  # Initialize self-loops
  if (add_self_loops) {
    diag(adj_mat) <- loop_value
  }

  # Normalize (sum by col must be 1)
  adj_mat_norm <- t(adj_mat / colSums(adj_mat))

  cur_mat_norm <- adj_mat_norm
  for (i_mcl in seq_len(max_iter)) {
    # message(i_mcl)
    last_mat <- cur_mat_norm

    exp_mat <- cur_mat_norm %^% mcl_expansion
    inf_mat <- exp_mat^mcl_inflation
    # inf_mat_norm <- t(inf_mat/colSums(inf_mat))
    inf_mat_norm <- apply(inf_mat, MARGIN = 2, FUN = function(matcol) {
      matcol / sum(matcol)
    })

    ## TODO: optional pruning?
    # idea: inspect matrix and set small values directly to zero (assume they would have reached there eventually anyways).
    if (identical(inf_mat_norm, last_mat)) {
      converged <- TRUE
      break
    }

    cur_mat_norm <- inf_mat_norm
  }

  if (converged & is.na(cur_mat_norm[1, 1])) {
    stop("An error occurred after convergence - maybe you set `add_self_loops` to FALSE?")
  }

  # getting the attractors - non-zero elements of the matrix diagonal
  clu_attractors <- which(diag(cur_mat_norm) > 0)
  # store the clusters
  clusters <- vector(mode = "list", length(clu_attractors))

  # the nodes in the same row as each attractor form a cluster
  for (i in seq_along(clu_attractors)) {
    cur_att <- clu_attractors[i]
    cur_members <- which(cur_mat_norm[cur_att, ] > 0)
    clusters[[i]] <- cur_members
  }

  # chop off the identical ones
  clusters <- unique(clusters)

  # from clusters sets to membership as label
  clu_memberships <- rep(NA, nrow(adj_mat))
  for (i in seq_along(clusters)) {
    this_cluster <- clusters[[i]]
    clu_memberships[this_cluster] <- i
  }

  # handling the singletons
  if (!allow_singletons) {
    dub <- duplicated(clu_memberships) + duplicated(clu_memberships, fromLast = TRUE)
    n_singletons <- sum(dub == 0)
    clu_memberships[dub == 0] <- 0
    # reshifting to not having gaps in the cluster numbers
    clu_memberships[dub != 0] <- clu_memberships[dub != 0] - n_singletons
  }

  res <- list()
  if (return_node_names & igraph::is_named(g)) {
    res$names <- V(g)$name
  }
  res$vcount <- vcount(g)
  res$algorithm <- "mcl"
  res$iterations <- i_mcl - 1
  res$membership <- clu_memberships
  if (return_esm) {
    res$esm <- cur_mat_norm
  }

  class(res) <- "communities" # to take advantage of the goodies for printing and co

  return(res)
}
