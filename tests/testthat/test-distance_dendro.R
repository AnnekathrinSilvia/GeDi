test_that("distance_dendro-no scores", {
  scores <- NULL
  expect_error(distance_dendro(scores))
  expect_error(distance_dendro(list()))
})

test_that("distance_dendro works", {
  distance_scores <- Matrix::Matrix(0.5, 20, 20)
  distance_scores[c(11:15), c(2:6)] <- 0.2
  dendro <- distance_dendro(distance_scores, cluster_method = "single")
  expect_type(dendro, "list")
})
