test_that("No genesets - getMeetMinMatrix", {
  genes <- list()
  expect_true(is.null(getMeetMinMatrix(genes, n_cores = 1)))
})

test_that("One empty geneset - getMeetMinMatrix", {
  genes <- list(list("PDHB", "VARS2"), list())
  m <- getMeetMinMatrix(genes, n_cores = 1)
  expect_equal(m[1, 2], 1)
})

test_that("Scoring identical sets - getMeetMinMatrix", {
  genes <- list(list("PDHB", "VARS2"), list("PDHB", "VARS2"))
  m <- getMeetMinMatrix(genes, n_cores = 1)
  expect_equal(m[1, 2], 0)
})

test_that("getMeetMinMatrix runs correctly", {
  genes <- list(list("PDHB", "VARS2"), list("IARS2", "PDHA1"))
  m <- getMeetMinMatrix(genes, n_cores = 1)
  expect_gte(m[1, 2], 0)
})
