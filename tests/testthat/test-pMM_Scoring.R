test_that("Empty Genesets - 1", {
  ppi <- Matrix::Matrix(0.5, 1, 1)
  expect_equal(sumInteraction(x = c(), y = c(), ppi = ppi), 0)
})

test_that("One empty Geneset - 1", {
  a <- c(1, 2)
  ppi <- Matrix::Matrix(0.5, 2, 2)
  expect_equal(sumInteraction(x = c(), y = a, ppi = ppi), 0)
  expect_equal(sumInteraction(x = a, y = c(), ppi = ppi), 0)
})

test_that("sumInteraction works", {
  ppi <- Matrix::Matrix(0.5, 2, 2)
  expect_equal(sumInteraction(x = c(1, 2), y = c(1), ppi = ppi), 1)
})

test_that("Empty Genesets - 2", {
  ppi <- Matrix::Matrix(0.5, 1, 1)
  maxInteract <- 0.5
  expect_equal(getInteractionScore(a = c(), ai = c(), b = c(), bi = c(), ppi = ppi, maxInteract = maxInteract), -1)
})

test_that("One empty Geneset - 2", {
  empty <- c()
  i <- c(1, 2)
  genes <- c("PDHB", "VARS2", "IARS2")
  ppi <- Matrix::Matrix(0.5, 2, 2)
  maxInteract <- 0.5
  expect_equal(getInteractionScore(a = empty, ai = i, b = genes, bi = i, ppi = ppi, maxInteract = maxInteract), -1)
  expect_equal(getInteractionScore(a = genes, ai = empty, b = genes, bi = i, ppi = ppi, maxInteract = maxInteract), -1)
  expect_equal(getInteractionScore(a = genes, ai = i, b = empty, bi = i, ppi = ppi, maxInteract = maxInteract), -1)
  expect_equal(getInteractionScore(a = genes, ai = i, b = genes, bi = empty, ppi = ppi, maxInteract = maxInteract), -1)
})

test_that("getInteractionScore works", {
  i <- c(1, 2, 3)
  genes <- c("PDHB", "VARS2", "IARS2")
  ppi <- Matrix::Matrix(0.5, 3, 3)
  maxInteract <- 0.5
  expect_gte(getInteractionScore(a = genes, ai = i, b = genes, bi = i, ppi = ppi, maxInteract = maxInteract), 0)
})

