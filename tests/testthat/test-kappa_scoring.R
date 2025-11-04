test_that("Empty genesets - getKappaMatrix", {
  expect_error(is.null(getKappaMatrix(genes = list())))
})

test_that("Scoring identical sets - getKappaMatrix", {
  genesets <- list(list("PHDB"), list("PHDB"))
  k <- getKappaMatrix(genesets)
  expect_equal(k[1, 2], 0)
})

test_that("getKappaMatrix runs correctly", {
  genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
  k <- getKappaMatrix(genesets)
  expect_equal(k[1, 2], 1)
})
