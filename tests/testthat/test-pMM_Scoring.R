test_that("Empty Genesets - getInteractionScore", {
  ppi <- data.frame(
    Gene1 = c("PDHB", "VARS2"),
    Gene2 = c("IARS2", "PDHA1"),
    combined_score = c(0.5, 0.2)
  )
  expect_equal(getInteractionScore(a = c(), b = c(), ppi = ppi, maxInteract = max(ppi$combined_score)), 0)
})

test_that("One empty Geneset - getInteractionScore", {
  a <- c("PDHB", "VARS2")
  ppi <- data.frame(
    Gene1 = c("PDHB", "VARS2"),
    Gene2 = c("IARS2", "PDHA1"),
    combined_score = c(0.5, 0.2)
  )
  maxInteract <- max(ppi$combined_score)
  expect_equal(getInteractionScore(a = c(), b = a, ppi = ppi, maxInteract = maxInteract), 0)
  expect_equal(getInteractionScore(a = a, b = c(), ppi = ppi, maxInteract = maxInteract), 0)
})

test_that("getInteractionScore runs correctly", {
  a <- c("PDHB", "VARS2")
  b <- c("IARS2", "PDHA1")
  ppi <- data.frame(
    Gene1 = c("PDHB", "VARS2"),
    Gene2 = c("IARS2", "PDHA1"),
    combined_score = c(0.5, 0.2)
  )
  maxInteract <- max(ppi$combined_score)
  expect_gte(getInteractionScore(a = a, b = b, ppi = ppi, maxInteract = maxInteract), 0)
})

test_that("Empty genesets - pMMlocal", {
  ppi <- data.frame()
  expect_equal(pMMlocal(a = c(), b = c(), ppi = ppi, maxInteract = 0), 1)
})

test_that("One geneset empty - pMMlocal", {
  a <- a <- c("PDHB", "VARS2")
  ppi <- data.frame()
  expect_equal(pMMlocal(a = a, b = c(), ppi = ppi, maxInteract = 0), 1)
})

test_that("pMMlocal runs correctly", {
  a <- c("PDHB", "VARS2")
  b <- c("IARS2", "PDHA1")
  ppi <- data.frame(
    Gene1 = c("PDHB", "VARS2"),
    Gene2 = c("IARS2", "PDHA1"),
    combined_score = c(0.5, 0.2)
  )
  maxInteract <- max(ppi$combined_score)
  expect_gte(pMMlocal(a = a, b = b, ppi = ppi, maxInteract = maxInteract), 0)
})

test_that("Empty genesets - getpMMMatrix", {
  ppi <- data.frame()
  expect_true(is.null(getpMMMatrix(genes = list(), ppi = ppi, n_cores = 1)))
})

test_that("getpMMMatrix runs correctly", {
  genes <- list(c("PDHB", "VARS2"), c("IARS2", "PDHA1"))
  ppi <- data.frame(
    Gene1 = c("PDHB", "VARS2"),
    Gene2 = c("IARS2", "PDHA1"),
    combined_score = c(0.5, 0.2)
  )
  m <- getpMMMatrix(genes = genes, ppi = ppi, n_cores = 1)
  expect_gte(m[1, 2], 0)
})
