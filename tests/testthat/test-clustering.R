test_that("checkInclusion - empty seeds", {
  seeds <- checkInclusion(list())
  expect_true(length(seeds) == 0)
  expect_identical(seeds, list())
})

test_that("checkInclusion - only one set", {
  seeds <- checkInclusion(list(c(1:20)))
  expect_identical(seeds, list(c(1:20)))

  seeds <- checkInclusion(list(c(1:20), c(1:20), c(1:20)))
  expect_identical(seeds, list(c(1:20)))
  expect_true(length(seeds) == 1)
})

test_that("checkInclusion runs correctly", {
  initial_seeds <- list(c(1:5), c(1:3), c(6:7), c(6:10))
  seeds <- checkInclusion(initial_seeds)
  expect_true(length(seeds) < length(initial_seeds))
  expect_identical(seeds, list(c(1:5), c(6:10)))
  expect_true(length(seeds) == 2)
})

test_that("seedFinding - no distance scores", {
  expect_true(is.null(seedFinding(list(), 0.3, 0.5)))
  expect_true(is.null(seedFinding(NULL, 0.3, 0.5)))
})

test_that("seedFinding - no distances smaller simthreshold", {
  distances <- Matrix::Matrix(0.9, 10, 10)
  seeds <- seedFinding(distances, simThreshold = 0.3, memThreshold = 0.5)
  expect_true(length(seeds) == 0)
  expect_type(seeds, "list")
})

test_that("seedFinding runs correctly", {
  distances <- Matrix::Matrix(0.9, 10, 10)
  seeds <- seedFinding(distances, simThreshold = 0.9, memThreshold = 0.5)
  expect_true(length(seeds) > 0)
  expect_type(seeds, "list")
})

test_that("clustering - no seeds", {
  expect_true(is.null(clustering(NULL, 0.5)))
  expect_true(length(clustering(list(), 0.5)) == 0)
})

test_that("clustering runs correctly", {
  seeds <- list(c(1:2), c(6:10), c(4:7))
  cluster <- clustering(seeds, 0.5)
  expect_true(length(cluster) <= length(seeds))
  expect_true(length(cluster) > 0)
  expect_type(cluster, "list")
})

test_that("getClusterDataTable - no geneset names", {
  expect_error(getClusterDatatable(list(), list()))
})

test_that("getClusterDatatable runs correctly", {
  geneset_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
  df <- getClusterDatatable(list(), geneset_names)
  expect_type(df, "list")
  expect_true(nrow(df) == 9)
  cluster <- list(c(1:4), c(5:8))
  expect_type(getClusterDatatable(cluster, geneset_names), "list")
})
