test_that("Empty genesets - getJaccardMatrix", {
  expect_error(is.null(getJaccardMatrix(genes = list())))
})

test_that("Scoring identical sets - getJaccardMatrix", {
  genesets <- list(list("PHDB"), list("PHDB"))
  k <- getJaccardMatrix(genesets)
  expect_equal(k[1, 2], 0)
})

test_that("getJaccardMatrix runs correctly", {
  genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
  k <- getJaccardMatrix(genesets)
  expect_equal(k[1, 2], 1)
})
