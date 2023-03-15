data(macrophage_topGO_example_small, package = "GeDi")
go_ids <- macrophage_topGO_example_small$Genesets
data(scores_macrophage_topGO_example_small, package = "GeDi")

test_that("No genesets - goSimilarity", {
  genes <- list()
  expect_equal(goSimilarity(geneset_ids = genes), -1)
})

test_that("Similarity calculation - no genesets", {
  expect_equal(goSimilarity(list()), -1)
  expect_error(goSimilarity(go_ids, method = "test"))
  expect_error(goSimilarity(go_ids, ontology = "test"))
  expect_error(goSimilarity(go_ids, species = "org.hs.eg.db"))
})

test_that("Similarity calculation runs correctly", {
  sim <- goSimilarity(geneset_ids = go_ids)
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
