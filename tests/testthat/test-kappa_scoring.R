test_that("no genesets", {
  genesets <- list()
  gene_names <- list()
  expect_equal(getKappaMatrix(genesets, gene_names), -1)
})

test_that("Scoring two sets", {
  genesets <- list(c("PHDB", "SOS2"), c("IARS2", "PDHA2"))
  gene_names <- list("A", "B")
  k <- getKappaMatrix(genesets, gene_names)
  expect_equal(k[1, 2], 1)
})

test_that("Scoring identical sets", {
  genesets <- list(c("PHDB"), c("PHDB"))
  gene_names <- list("A", "B")
  k <- getKappaMatrix(genesets, gene_names)
  expect_equal(k[1, 2], 0)
})
