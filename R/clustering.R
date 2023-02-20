#' Check for subset inclusion
#'
#' Check if every set is unique and there are no subset included in larger sets
#'
#' @param seeds A `list` of sets of numerical values
#'
#' @return A `list` of unique sets
#' @export
#' @importFrom rje is.subset
#'
#' @examples
#' seeds <- list(c(1:5), c(2:5), c(6:10))
#' s <- checkInclusion(seeds)
checkInclusion <- function(seeds) {
  remove <- c()
  if (length(seeds) < 2) {
    return(seeds)
  }
  for (i in 1:(length(seeds) - 1)) {
    s1 <- seeds[[i]]
    for (j in (i + 1):length(seeds)) {
      s2 <- seeds[[j]]
      if (rje::is.subset(s1, s2)) {
        remove <- c(remove, i)
      } else if (rje::is.subset(s2, s1)) {
        remove <- c(remove, j)
      }
    }
  }
  remove <- unique(remove)
  if (length(remove) == 0) {
    return(seeds)
  } else{
    return(seeds[-remove])
  }
}

#' Find clustering seeds
#'
#' Find initial seeds for the clustering. #TODO: Verlinken auf die David publication
#'
#' @param distances A [Matrix::Matrix()] of (distance) scores
#' @param simThreshold Numerical value, a threshold of what is considered a
#'                     close relationsship between two elements in distances
#' @param memThreshold Numerical value, a threshold used to identify members of
#'                     a seeds
#'
#' @return A `list` of seeds which can be used for clustering
#' @export
#'
#' @examples
#' m <- Matrix::Matrix(runif(100, min = 0, max = 1), 10, 10)
#' seeds <- seedFinding(distances = m, simThreshold = 0.3, memThreshold = 0.5)
seedFinding <- function(distances, simThreshold, memThreshold) {
  if (is.null(distances) || length(distances) == 0) {
    return(NULL)
  }
  seeds <- list()

  reach <-
    apply(distances, 1, function(x)
      as.numeric(x <= simThreshold))
  for (i in 1:nrow(distances)) {
    if (sum(reach[i, ], na.rm = TRUE) >= 2) {
      members <- which(reach[i,] == 1)
      includethreshold <-
        (length(members) ^ 2 - length(members)) * memThreshold
      reach_red <- reach[members, members]
      in_reach <- sum(reach_red)
      if (in_reach >= includethreshold) {
        members <- c(members, i)
        seeds <- c(list(sort(members)), seeds)
      }
    }
  }
  seeds <- checkInclusion(unique(seeds))
  seeds <- lapply(seeds, unique)
  return(seeds)
}


#' Find cluster from initial seeds
#'
#' Find clusters from initial seeds according to (#TODO: David publication linken)
#'
#' @param seeds A `list` of seeds, e.g. determined by seedFinding
#' @param threshold Numerical value, threshold for merging seeds
#'
#' @return A `list` of clusters
#' @export
#'
#' @examples
#' seeds <- list(c(1:5), c(6:10))
#' cluster <- clustering(seeds, 0.5)
clustering <- function(seeds, threshold) {
  if (length(seeds) <= 1) {
    return(seeds)
  }
  mergeable <- rep(TRUE, length(seeds))
  while (any(mergeable)) {
    index <- which(mergeable)[1]
    if (index > length(seeds)) {
      break
    }
    s1 <- seeds[[index]]
    l <- length(seeds)
    for (j in 1:length(seeds)) {
      s2 <- seeds[[j]]
      int <- intersect(s1, s2)
      union <- union(s1, s2)
      if (length(int) >= (threshold * length(union))) {
        remove <- list(s1, s2)
        seeds <- seeds[!(seeds %in% remove)]
        seeds <- c(list(union), seeds)
        mergeable <- mergeable[-c(index, j)]
        mergeable <- c(TRUE, mergeable)
        break
      }
    }
    if (l == length(seeds)) {
      mergeable[[index]] <- FALSE
    }

  }
  return(seeds)
}

#' Get Dataframe of Cluster
#'
#' Get Dataframe from a `list` of clusters
#'
#' @param cluster A `list` of clusters
#' @param geneset_names A vector of geneset names
#'
#' @return A `data.frame` mapping each geneset to the cluster(s) it belongs to
#' @export
#'
#' @examples
#' cluster <- list(c(1:5), c(4:8))
#' geneset_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
#' df <- getClusterDatatable(cluster, geneset_names)
getClusterDatatable <- function(cluster, geneset_names) {
  stopifnot(length(geneset_names) > 0)
  n_gs <- length(geneset_names)
  df <- vector("list", n_gs)

  if (length(cluster) == 0) {
    df <- data.frame(Cluster = rep("No associated Cluster", n_gs))
    rownames(df) <- geneset_names
    return(df)
  }

  for (i in 1:length(cluster)) {
    for (j in cluster[[i]]) {
      entry <- df[[j]]
      if (is.null(entry)) {
        entry <- i
      } else{
        entry <- c(entry, i)
      }

      df[[j]] <- entry
    }
  }

  df <- data.frame(matrix(df, nrow = n_gs, ncol = 1))
  colnames(df) <- c("Cluster")
  cluster <- df$Cluster
  cluster <-
    lapply(cluster, function(x)
      if (is.null(x)) {
        "No associated Cluster"
      } else{
        x
      })
  df$Cluster <- cluster
  rownames(df) <- geneset_names
  return(df)
}
