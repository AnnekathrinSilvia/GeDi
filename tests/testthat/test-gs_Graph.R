data(scores_macrophage_topGO_example_small, package = "GeDi")
cluster_macrophage_topGO_example_small <- list(c(1, 3, 9, 41), c(15, 28))
data(macrophage_topGO_example_small, package = "GeDi")
gs_names <- macrophage_topGO_example_small$Genesets
data(macrophage_Reactome_example, package = "GeDi")
data(macrophage_KEGG_example, package = "GeDi")
gs_names_reactome <- macrophage_Reactome_example$ID
gs_names_reactome <- gs_names_reactome[1:nrow(scores_macrophage_topGO_example_small)]
gs_names_kegg <- macrophage_KEGG_example$ID
gs_names_kegg <- gs_names_kegg[1:nrow(scores_macrophage_topGO_example_small)]
genes <- prepareGenesetData(macrophage_topGO_example_small)

test_that("getAdjacencyMatrix - No distance scores", {
  expect_true(is.null(getAdjacencyMatrix(NULL, 0.5)))
})

test_that("getAdjacencyMatrix runs correctly", {
  m <- getAdjacencyMatrix(scores_macrophage_topGO_example_small, 0.3)
  expect_s4_class(m, "Matrix")
  expect_equal(max(m), 1)
  expect_equal(min(m), 0)
})

test_that("buildGraph runs correctly", {
  m <- getAdjacencyMatrix(scores_macrophage_topGO_example_small, 0.3)
  g <- buildGraph(m)
  expect_type(g, "list")

  scores_reactome <- scores_macrophage_topGO_example_small
  rownames(scores_reactome) <- colnames(scores_reactome) <- gs_names_reactome
  g <- buildGraph(getAdjacencyMatrix(scores_reactome, 0.3))
  expect_type(g, "list")

  scores_kegg <- scores_macrophage_topGO_example_small
  rownames(scores_kegg) <- colnames(scores_kegg) <- gs_names_kegg
  g <- buildGraph(getAdjacencyMatrix(scores_kegg, 0.3))
  expect_type(g, "list")
})

test_that("getClusterAdjacencyMatrix runs correctly", {
  adj <- getClusterAdjacencyMatrix(
    cluster_macrophage_topGO_example_small,
    gs_names
  )
  expect_equal(max(adj), 1)
  expect_equal(min(adj), 0)
  expect_s4_class(adj, "Matrix")
})

test_that("getClusterAdjacencyMatrix runs with no cluster", {
  cluster <- list()
  adj <- getClusterAdjacencyMatrix(cluster, gs_names)
  expect_equal(max(adj), 0)
  expect_equal(min(adj), 0)
  expect_s4_class(adj, "Matrix")
})

test_that("getClusterAdjacencyMatrix - no gs_names", {
  expect_error(getClusterAdjacencyMatrix(
    cluster_macrophage_topGO_example_small,
    list()
  ))
  expect_error(getClusterAdjacencyMatrix(
    cluster_macrophage_topGO_example_small,
    list("a", "b")
  ))
})

test_that("getBipartiteGraph - one input missing", {
  expect_error(getBipartiteGraph(
    cluster = list(),
    geneset_names = list(),
    genes = list()
  ))

  expect_error(getBipartiteGraph(
    cluster_macrophage_topGO_example_small,
    list(),
    genes
  ))

  expect_error(getBipartiteGraph(
    cluster_macrophage_topGO_example_small,
    gs_names,
    list()
  ))
})

test_that("getBipartiteGraph runs correctly", {
  g <- getBipartiteGraph(
    cluster_macrophage_topGO_example_small,
    gs_names,
    genes
  )
  expect_type(g, "list")
})

test_that("buildClusterGraph - no cluster", {
  cluster <- list()

  expect_error(  graph <- buildClusterGraph(
    cluster,
    macrophage_topGO_example_small,
    gs_names
  ))
})

test_that("buildClusterGraph coloring", {
  expect_error(buildClusterGraph(
    cluster_macrophage_topGO_example_small,
    macrophage_topGO_example_small,
    gs_names,
    color_by = "test_column"
  ))

  graph <- buildClusterGraph(
    cluster_macrophage_topGO_example_small,
    macrophage_topGO_example_small,
    gs_names,
    color_by = "Annotated"
  )

  graph <- buildClusterGraph(
    cluster_macrophage_topGO_example_small,
    macrophage_topGO_example_small,
    gs_names,
    color_by = "p.value_elim"
  )
  graph <- buildClusterGraph(
    cluster_macrophage_topGO_example_small,
    macrophage_topGO_example_small,
    gs_names,
    color_by = "Significant"
  )
  expect_type(graph, "list")
})

test_that("buildClusterGraph works correctly", {
  graph <- buildClusterGraph(
    cluster_macrophage_topGO_example_small,
    macrophage_topGO_example_small,
    gs_names
  )
  expect_type(graph, "list")
})

test_that(".graphMetricsGenesetsDT runs correctly", {
  m <- getAdjacencyMatrix(scores_macrophage_topGO_example_small, 0.3)
  g <- buildGraph(m)
  dt <- .graphMetricsGenesetsDT(
    g,
    macrophage_topGO_example_small
  )
  expect_type(dt, "list")
})
