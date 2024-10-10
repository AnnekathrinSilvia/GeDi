data(macrophage_topGO_example_small, package = "GeDi")
go_ids <- macrophage_topGO_example_small$Genesets
data(scores_macrophage_topGO_example_small, package = "GeDi")

test_that("No genesets - goDistance", {
  genes <- list()
  expect_equal(goDistance(geneset_ids = genes), -1)
})

test_that("Similarity calculation - no genesets", {
  expect_equal(goDistance(list()), -1)
  expect_error(goDistance(go_ids, method = "test"))
  expect_error(goDistance(go_ids, ontology = "test"))
  expect_error(goDistance(go_ids, species = "org.hs.eg.db"))
})

test_that("Similarity calculation runs correctly", {
  sim <- goDistance(geneset_ids = go_ids)
  expect_gte(sim[1, 1], 0)
})

test_that("Scaling runs correctly", {
  expect_error(scaleGO(
    scores = scores_macrophage_topGO_example_small,
    geneset_ids = list()
  ))
  scaled <- scaleGO(
    scores = scores_macrophage_topGO_example_small,
    geneset_ids = go_ids
  )
  expect_gte(scaled[1, 1], 0)
})
