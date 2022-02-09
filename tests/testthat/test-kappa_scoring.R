test_that("no genesets", {
  genesets <- list()
  expect_equal(getKappaMatrix(genesets), -1)
})

test_that("Correct Scoring", {
  genesets <- list(c("PDHB", "VARS2", "IARS2"), c("IARS2", "PDHA2"))
  k <- getKappaMatrix(genesets)
  expect_gte(k[1, 1], 0)
})
