test_that("Empty genesets - calculateSorensenDice", {
  expect_equal(calculateSorensenDice(a = c(), b = c()), 1)
})

test_that("One empty geneset - calculateSorensenDice", {
  a <- c("PDHB", "VARS2")
  b <- c("IARS2", "PDHA1")
  expect_equal(calculateSorensenDice(a = a, b = c()), 1)
  expect_equal(calculateSorensenDice(a = c(), b = b), 1)
})

test_that("calculateSorensenDice runs correctly", {
  a <- c("PDHB", "VARS2")
  b <- c("IARS2", "PDHA1")
  expect_gte(calculateSorensenDice(a = a, b = b), 0)
})

test_that("Empty genesets - getSorensenDiceMatrix", {
  expect_true(is.null(getSorensenDiceMatrix(genes = list(), n_cores = 1)))
})

test_that("Scoring identical sets - getSorensenDiceMatrix", {
  genesets <- list(list("PHDB"), list("PHDB"))
  k <- getSorensenDiceMatrix(genesets, n_cores = 1)
  expect_equal(k[1, 2], 0)
})

test_that("getSorensenDiceMatrix runs correctly", {
  genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
  k <- getSorensenDiceMatrix(genesets, n_cores = 1)
  expect_equal(k[1, 2], 1)
})
