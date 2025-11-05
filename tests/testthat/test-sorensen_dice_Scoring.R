test_that("Empty genesets - getSorensenDiceMatrix", {
  expect_error(getSorensenDiceMatrix(genes = list()))
})

test_that("Scoring identical sets - getSorensenDiceMatrix", {
  genesets <- list(list("PHDB"), list("PHDB"))
  k <- getSorensenDiceMatrix(genesets)
  expect_equal(k[1, 2], 0)
})

test_that("getSorensenDiceMatrix runs correctly", {
  genesets <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
  k <- getSorensenDiceMatrix(genesets)
  expect_equal(k[1, 2], 1)
})
