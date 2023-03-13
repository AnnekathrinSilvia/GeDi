test_that("getAdjacencyMatrix - No distance scores", {
  expect_true(is.null(getAdjacencyMatrix(NULL, 0.5)))
})

test_that("getAdjacencyMatrix runs correctly", {
  distance_scores <- Matrix::Matrix(stats::runif(1000, min = 0, max = 1), 100, 100)
  threshold <- 0.3
  m <- getAdjacencyMatrix(distance_scores, threshold)
  expect_s4_class(m, "Matrix")
  expect_equal(max(m), 1)
  expect_equal(min(m), 0)
})

test_that("buildGraph runs correctly", {
  distance_scores <- Matrix::Matrix(stats::runif(1000, min = 0, max = 1), 100, 100)
  threshold <- 0.3
  geneset_names <- as.character(stats::runif(100, min = 0, max = 1))
  rownames(distance_scores) <- colnames(distance_scores) <- geneset_names
  m <- getAdjacencyMatrix(distance_scores, threshold)
  g <- buildGraph(m)
  expect_type(g, "list")
})

test_that("getClusterAdjacencyMatrix runs correctly", {
  cluster <- list(c(10:25))
  geneset_names <- as.character(stats::runif(29, min = 0, max = 1))
  adj <- getClusterAdjacencyMatrix(cluster, geneset_names)
  expect_equal(max(adj), 1)
  expect_equal(min(adj), 0)
  expect_s4_class(adj, "Matrix")
})

test_that("getClusterAdjacencyMatrix runs with no cluster", {
  cluster <- list()
  geneset_names <- as.character(stats::runif(29, min = 0, max = 1))
  adj <- getClusterAdjacencyMatrix(cluster, geneset_names)
  expect_equal(max(adj), 0)
  expect_equal(min(adj), 0)
  expect_s4_class(adj, "Matrix")
})

test_that("getClusterAdjacencyMatrix - no gs_names", {
  cluster <- list(c(1:5), c(6:9, 1))
  gs_names <- list("a", "b")
  expect_error(getClusterAdjacencyMatrix(cluster, list()))
  expect_error(getClusterAdjacencyMatrix(cluster, gs_names))
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
  geneset_names <- as.character(stats::runif(9, min = 0, max = 1))
  expect_error(getBipartiteGraph(cluster, geneset_names, list()))
})

test_that("getBipartiteGraph runs correctly", {
  cluster <- list(c(1:2))
  geneset_names <- c("a", "b")
  genes <- list(c("PDHB", "VARS2"), c("IARS2", "PDHA1"))
  g <- getBipartiteGraph(cluster, geneset_names, genes)
  expect_type(g, "list")
})

test_that("buildClusterGraph - no cluster", {
  cluster <- list()
  genes <- list(
    c("PDHB", "VARS2"), c("IARS2", "PDHA1"),
    c("AAAS", "ABCE1"), c("ABI1", "AAR2"), c("AATF", "AMFR"),
    c("BMS1", "DAP3"), c("AURKAIP1", "CHCHD1"), c("IARS2"),
    c("AHI1", "ALMS1")
  )
  gs_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
  geneset_df <- data.frame(
    Genesets = gs_names
  )
  geneset_df$Genes <- genes
  graph <- buildClusterGraph(cluster,
                             geneset_df,
                             genes)
  expect_type(graph, "list")
  expect_equal(gorder(graph), 0)
  expect_equal(gsize(graph), 0)
})

test_that("buildClusterGraph works correctly", {
  cluster <- list(c(1:5), c(6:9, 1))
  genes <- list(
    c("PDHB", "VARS2"), c("IARS2", "PDHA1"),
    c("AAAS", "ABCE1"), c("ABI1", "AAR2"), c("AATF", "AMFR"),
    c("BMS1", "DAP3"), c("AURKAIP1", "CHCHD1"), c("IARS2"),
    c("AHI1", "ALMS1")
  )
  gs_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
  geneset_df <- data.frame(
    Genesets = gs_names
  )
  geneset_df$Genes <- genes
  graph <- buildClusterGraph(cluster,
                             geneset_df,
                             genes)
  expect_type(graph, "list")
})
