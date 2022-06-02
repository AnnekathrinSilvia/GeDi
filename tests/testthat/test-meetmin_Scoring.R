test_that("no genesets", {
  genes <- list()
  genes_names <- list()
  expect_equal(getMeetMinMatrix(genes, genes_names), -1)
})

test_that("Scoring identical sets", {
  genes <- list(c("IARS2"), c("IARS2"))
  genes_names <- c("A", "B")
  m <- getMeetMinMatrix(genes, genes_names)
  expect_equal(m[1, 1], 0)
  })

test_that("Scoring two sets", {
  genes <- list(c("PDHB", "SOS2"), c("IARS2", "LIPG"))
  genes_names <- c("A", "B")
  m <- getMeetMinMatrix(genes, genes_names)
  expect_equal(m[1, 2], 1)
})
