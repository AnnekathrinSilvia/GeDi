#' Title
#'
#' @param seeds
#'
#' @return
#' @export
#'
#' @examples
checkInclusion <- function(seeds) {
  remove <- c()
  if(length(seeds) < 2){
    return(seeds)
  }
  for (i in 1:(length(seeds) - 1)) {
    s1 <- seeds[[i]]
    l1 <- length(s1)
    for (j in (i + 1):length(seeds)) {
      s2 <- seeds[[j]]
      l2 <- length(s2)
      if (l1 < l2) {
        if (length(intersect(s1, s2)) == l1) {
          remove <- c(remove, i)
          break
        }
      } else{
        if (length(intersect(s2, s1)) == l2) {
          remove <- c(remove, j)
        }
      }
    }
  }
  remove <- unique(remove)
  if(length(remove) == 0){
    return(seeds)
  }else{
    return(seeds[-remove])
  }
}

#' Title
#'
#' @param distances
#' @param simThreshold
#' @param memThreshold
#'
#' @return
#' @export
#'
#' @examples
seedFinding <- function(distances, simThreshold, memThreshold){
  # simthreshold: what is considered 'close' relationship
  # memthreshold: how many members of a possible seed need a close relationship for the seed to be considered
  seeds <- list()

  # Build matrix with xij = 1 indicating i and j are close (sim(i, j) <= simThreshold)
  reach <- apply(distances, 1, function(x) as.numeric(x <= simThreshold))
  for(i in 1:nrow(distances)){
    if(sum(reach[i,], na.rm = TRUE) >= 2){
      members <- which(reach[i, ] == 1)
      includethreshold <- (length(members)^2 - length(members)) * memThreshold
      reach_red <- reach[members, members]
      in_reach <- sum(reach_red)
      if(in_reach >= includethreshold){
        members <- c(members, i)
        seeds <- c(list(sort(members)), seeds)
      }
    }
  }
  seeds <- checkInclusion(unique(seeds))
  seeds <- lapply(seeds, unique)
  return(seeds)
}


#' Title
#'
#' @param seeds
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
clustering <- function(seeds, threshold){
  if(length(seeds) <= 1){
    return(seeds)
  }
  mergeable <- rep(TRUE, length(seeds))
  while(any(mergeable)){
    index <- which(mergeable)[1]
    if(index > length(seeds)){
      break
    }
    s1 <- seeds[[index]]
    l <- length(seeds)
    for(j in 1:length(seeds)){
      s2 <- seeds[[j]]
      int <- intersect(s1, s2)
      union <- union(s1, s2)
      if(length(int) >= (threshold * length(union))){
        remove <- list(s1, s2)
        seeds <- seeds[!(seeds %in% remove)]
        seeds <- c(list(union), seeds)
        mergeable <- mergeable[-c(index, j)]
        mergeable <- c(TRUE, mergeable)
        break
      }
    }
    if(l == length(seeds)){
      mergeable[[index]] <- FALSE
    }

  }
  return(seeds)
}

#' Title
#'
#' @param cluster
#' @param n_gs
#'
#' @return
#' @export
#'
#' @examples
getClusterDatatable <- function(cluster, geneset_names){
  n_gs <- length(geneset_names)
  df <- vector("list", n_gs)

  if(length(cluster) == 0){
    df <- data.frame(Cluster = rep("No associated Cluster", n_gs))
    rownames(df) <- geneset_names
    return(df)
  }

  for(i in 1:length(cluster)){
    for(j in cluster[[i]]){
      entry <- df[[j]]
      if(is.null(entry)){
        entry <- i
      }else{
        entry <- c(entry, i)
      }

      df[[j]] <- entry
    }
  }

  df <- data.frame(matrix(df, nrow = n_gs, ncol = 1))
  colnames(df) <- c("Cluster")
  cluster <- df$Cluster
  cluster <- lapply(cluster, function(x) if(is.null(x)){"No associated Cluster"}else{x})
  df$Cluster <- cluster
  rownames(df) <- geneset_names
  return(df)
}


