test_that("getAdjacencyMatrix - No distance scores", {
  expect_true(is.null(getAdjacencyMatrix(NULL, 0.5)))
})

test_that("getAdjacencyMatrix runs correctly", {
  distance_scores <- Matrix::Matrix(runif(1000, min = 0, max = 1), 100, 100)
  threshold <- 0.3
  m <- getAdjacencyMatrix(distance_scores, threshold)
  expect_s4_class(m, "Matrix")
  expect_equal(max(m), 1)
  expect_equal(min(m), 0)
})

test_that("buildGraph runs correctly", {
  distance_scores <- Matrix::Matrix(runif(1000, min = 0, max = 1), 100, 100)
  threshold <- 0.3
  geneset_names <- as.character(runif(100, min = 0, max = 1))
  rownames(distance_scores) <- colnames(distance_scores) <- geneset_names
  m <- getAdjacencyMatrix(distance_scores, threshold)
  g <- buildGraph(m)
  expect_type(g, "list")
})

test_that("getClusterAdjacencyMatrix runs correctly", {
  cluster <- list(c(10:25))
  geneset_names <- as.character(runif(29, min = 0, max = 1))
  adj <- getClusterAdjacencyMatrix(cluster, geneset_names)
  expect_equal(max(adj), 1)
  expect_equal(min(adj), 0)
  expect_s4_class(adj, "Matrix")
})

test_that("getClusterAdjacencyMatrix runs with no cluster", {
  cluster <- list()
  geneset_names <- as.character(runif(29, min = 0, max = 1))
  adj <- getClusterAdjacencyMatrix(cluster, geneset_names)
  expect_equal(max(adj), 0)
  expect_equal(min(adj), 0)
  expect_s4_class(adj, "Matrix")
})

test_that("getBipartiteGraph - no clusters", {
  expect_error(getBipartiteGraph(cluster = list(), geneset_names = list(), genes = list()))
})

test_that("getBipartiteGraph - no geneset names", {
  cluster <- list(c(1:20), c(4:6))
  genes <- list(c("PDHB", "VARS2"), c("IARS2", "PDHA1"))
  expect_error(getBipartiteGraph(cluster, list(), genes))
})

test_that("getBipartiteGraph - no genes", {
  cluster <- list(c(1:5), c(8:9))
  geneset_names <- as.character(runif(9, min = 0, max = 1))
  expect_error(getBipartiteGraph(cluster, geneset_names, list()))
})

test_that("getBipartiteGraph runs correctly", {
  cluster <- list(c(1:2))
  geneset_names <- c("a", "b")
  genes <- list(c("PDHB", "VARS2"), c("IARS2", "PDHA1"))
  g <- getBipartiteGraph(cluster, geneset_names, genes)
  expect_type(g, "list")
})
