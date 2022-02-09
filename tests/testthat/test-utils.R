test_that("no genesets", {
  genesets <- list()
  expect_equal(length(getGenes(genesets = genesets)), 0)
})

test_that("no genesets", {
  genes <- list()
  expect_equal(length(getGenesIndexed(genes = genes)), 0)
})
