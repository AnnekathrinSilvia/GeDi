distances <- Matrix::Matrix(0.9, 10, 10)
rownames(distances) <- colnames(distances) <- c("a", "b", "c", "d", "e",
                                                "f", "g", "h", "i", "j")

test_that("checkInclusion - empty seeds", {
  seeds <- checkInclusion(list())
  expect_true(length(seeds) == 0)
  expect_identical(seeds, list())

  seeds <- checkInclusion(list(c(), c()))
  expect_true(length(seeds) == 0)
  expect_identical(seeds, list())
})

test_that("checkInclusion - only one set", {
  seeds <- checkInclusion(list(c(1:20)))
  expect_identical(seeds, list(c(1:20)))

  seeds <- checkInclusion(list(c(1:20), c(1:20), c(1:20)))
  expect_identical(seeds, list(c(1:20)))
  expect_true(length(seeds) == 1)

  seeds <- checkInclusion(list(c(1:20), c()))
  expect_identical(seeds, list(c(1:20)))
  expect_true(length(seeds) == 1)

  seeds <- checkInclusion(list(c(1:20), c(20:1)))
  expect_identical(seeds, list(c(20:1)))
  expect_true(length(seeds) == 1)

})

test_that("checkInclusion - no seeds to remove", {
  seeds <- c(list(c(1:5)), list(6:10))
  expect_identical(seeds, checkInclusion(seeds))
})

test_that("checkInclusion runs correctly", {
  initial_seeds <- list(c(1:5), c(1:3), c(6:7), c(6:10))
  seeds <- checkInclusion(initial_seeds)
  expect_true(length(seeds) < length(initial_seeds))
  expect_identical(seeds, list(c(1:5), c(6:10)))
  expect_true(length(seeds) == 2)
})

test_that("checkInclusion works with Strings", {
  initial_seeds <- list(c("Gene1"), c("Gene2"), c("Gene1", "Gene3"), c("Gene2", "Gene4"))
  seeds <- checkInclusion(initial_seeds)
  expect_true(length(seeds) < length(initial_seeds))
  expect_identical(seeds, list(c("Gene1", "Gene3"), c("Gene2", "Gene4")))
  expect_true(length(seeds) == 2)

  # Check case-sensitivity
  initial_seeds <- list(c("gene1"), c("Gene2"), c("Gene1", "Gene3"), c("Gene2", "Gene4"))
  seeds <- checkInclusion(initial_seeds)
  expect_true(length(seeds) < length(initial_seeds))
  expect_identical(seeds, list(c("gene1"), c("Gene1", "Gene3"), c("Gene2", "Gene4")))
  expect_true(length(seeds) == 3)
})

test_that("seedFinding - no distance scores", {
  expect_true(is.null(seedFinding(list(), 0.3, 0.5)))
  expect_true(is.null(seedFinding(NULL, 0.3, 0.5)))

  # All distance scores are Nan
  distances <- Matrix::Matrix(NaN, 10, 10)
  expect_equal(seedFinding(distances, 0.3, 0.5), list())
})

test_that("seedFinding - no distances smaller simthreshold", {
  seeds <- seedFinding(distances, simThreshold = 0.3, memThreshold = 0.5)
  expect_true(length(seeds) == 0)
  expect_type(seeds, "list")
})

test_that("seedFinding runs correctly", {
  seeds <- seedFinding(distances, simThreshold = 0.9, memThreshold = 0.5)
  expect_true(length(seeds) > 0)
  expect_type(seeds, "list")

  # all genes in one seed
  seeds <- seedFinding(distances, simThreshold = 1, memThreshold = 0)
  expect_true(length(seeds) == 1)
  expect_type(seeds, "list")
})

test_that("fuzzyClustering - no seeds", {
  expect_true(is.null(fuzzyClustering(NULL, 0.5)))
  expect_true(length(fuzzyClustering(list(), 0.5)) == 0)
})

test_that("fuzzyClustering runs correctly", {
  # Only one seed
  seeds <- list(c(1:2))
  cluster <- fuzzyClustering(seeds, 0.5)
  expect_true(length(cluster) == length(seeds))
  expect_true(length(cluster) == 1)
  expect_type(cluster, "list")
  expect_equal(cluster, seeds)

  # several seeds
  seeds <- list(c(1:2), c(6:10), c(4:7))
  cluster <- fuzzyClustering(seeds, 0.5)
  expect_true(length(cluster) <= length(seeds))
  expect_true(length(cluster) > 0)
  expect_type(cluster, "list")

  # several seeds resulting in one cluster
  seeds <- list(c(1:6), c(3:6), c(4:6))
  cluster <- fuzzyClustering(seeds, 0.5)
  expect_true(length(cluster) <= length(seeds))
  expect_true(length(cluster) == 1)
  expect_type(cluster, "list")
})


test_that("clustering - wrong cluster_method", {
  expect_error(clustering(NULL, 0.3, "test"))
  expect_error(clustering(NULL, 0.3, "Markov"))
  expect_error(clustering(NULL, 0.3, "Louvain"))
})

test_that("clustering - no distance scores", {
  expect_error(clustering(NULL, 0.3))
  expect_error(clustering(list(), 0.3))
})

test_that("clustering works correctly", {
  # No scores under threshold, no clusters expected
  cluster <- clustering(distances, 0.3)
  expect_equal(cluster, list())
  expect_true(length(cluster) == 0)


  # Test cluster formation
  cluster <- clustering(distances, 1)
  expect_true(length(cluster) > 0)
  # Test cluster formation using Markov clustering
  cluster <- clustering(distances, 1, cluster_method = "markov")
  expect_true(length(cluster) > 0)
})

test_that("kNN clustering - no distance scores", {
  expect_true(is.null(kNN_clustering(NULL, 3)))
})

test_that("kNN clustering works correctly", {
  knn <- kNN_clustering(distances, 3)
  expect_true(length(knn) > 0 )
})

test_that("getClusterDataTable - no geneset names", {
  expect_error(.getClusterDatatable(list(), list()))
})

test_that("getClusterDatatable runs correctly", {
  geneset_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
  df <- .getClusterDatatable(list(), geneset_names)
  expect_type(df, "list")
  expect_true(nrow(df) == 9)
  cluster <- list(c(1:5), c(5:8))
  expect_type(.getClusterDatatable(cluster, geneset_names), "list")
})
