test_that("no genesets", {
  genes <- list()
  expect_equal(getMeetMinMatrix(genes), -1)
})

test_that("scoring works", {
  genes <- list(c("IARS2"), c("IARS2"))
  m <- getMeetMinMatrix(genes)
  expect_equal(m[1, 1], 0)
  })
