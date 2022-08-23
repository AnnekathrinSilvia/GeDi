test_that("No genesets - goSimilarity", {
  genes <- list()
  expect_equal(goSimilarity(genesets = genes), -1)
})

#test_that("Similarity calculation works",{
#  genesets <- list(c("PDHB", "VARS2", "IARS2"), c("IARS2", "PDHA2"))
#  sim <- goSimilarity(genesets = genesets)
#  expect_gte(sim[1, 1], 0)
#})

#test_that("Scaling works", {
#  genesets <- list(c("PDHB", "VARS2", "IARS2"), c("IARS2", "PDHA2"))
#  scores <- Matrix::Matrix(1, 1, 1)
#  scaled <- scaleGO(scores = scores, genesets = genesets)
#  expect_gte(scaled[1, 1], 0)
#})
