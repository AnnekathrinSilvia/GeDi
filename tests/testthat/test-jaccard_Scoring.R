test_that("Empty genesets - calculateJaccard", {
  expect_equal(calculateJaccard(a = c(), b = c()), 1)
})

test_that("One empty geneset - calculateJaccard", {
  a <- c("PDHB", "VARS2")
  b <- c("IARS2", "PDHA1")
  expect_equal(calculateJaccard(a = a, b = c()), 1)
  expect_equal(calculateJaccard(a = c(), b = b), 1)
})

test_that("calculateJaccard runs correctly", {
  a <- c("PDHB", "VARS2")
  b <- c("IARS2", "PDHA1")
  expect_gte(calculateJaccard(a = a, b = b), 0)
})

test_that("Empty genesets - getJaccardMatrix", {
  expect_true(is.null(getJaccardMatrix(genes = list(), n_cores = 1)))
})

test_that("Scoring identical sets - getJaccardMatrix", {
  genesets <- list(list("PHDB"), list("PHDB"))
  k <- getJaccardMatrix(genesets, n_cores = 1)
  expect_equal(k[1, 2], 0)
})

test_that("getJaccardMatrix runs correctly", {
  genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
  k <- getJaccardMatrix(genesets, n_cores = 1)
  expect_equal(k[1, 2], 1)
})
