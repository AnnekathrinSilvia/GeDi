test_that("Empty genesets - calculateKappa", {
  all_genes <- c("PDHB", "VARS2", "IARS2", "PDHA1")
  expect_equal(calculateKappa(a = c(), b = c(), all_genes = all_genes), 1)
})

test_that("One empty geneset - calculateKappa", {
  a <- c("PDHB", "VARS2")
  b <- c("IARS2", "PDHA1")
  all_genes <- c("PDHB", "VARS2", "IARS2", "PDHA1")
  expect_equal(calculateKappa(a = a, b = c(), all_genes = all_genes), 1)
  expect_equal(calculateKappa(a = c(), b = b, all_genes = all_genes), 1)
})

test_that("No unique genes - calculateKappa", {
  a <- c("PDHB", "VARS2")
  b <- c("IARS2", "PDHA1")
  expect_error(calculateKappa(a = c(), b = c(), all_genes = c()))
  expect_error(calculateKappa(a = a, b = c(), all_genes = c()))
  expect_error(calculateKappa(a = c(), b = b, all_genes = c()))
  expect_error(calculateKappa(a = a, b = b, all_genes = c()))
})

test_that("calculateKappa runs correctly", {
  a <- c("PDHB", "VARS2")
  b <- c("IARS2", "PDHA1")
  all_genes <- c("PDHB", "VARS2", "IARS2", "PDHA1")
  expect_gte(abs(calculateKappa(a = a, b = b, all_genes = all_genes)), 0)
})

test_that("Empty genesets - getKappaMatrix", {
  expect_true(is.null(getKappaMatrix(genes = list())))
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
