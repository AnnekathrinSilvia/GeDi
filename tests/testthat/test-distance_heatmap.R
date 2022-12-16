test_that("No distance scores", {
  distance_scores <- NULL
  expect_error(distance_heatmap(distance_scores = distance_scores))
})

test_that("Heatmap is correctly created", {
  distance_scores <- Matrix::Matrix(runif(1000, min=0, max=1), 100, 100)
  p <- distance_heatmap(distance_scores)
  expect_type(p, "list")
})

test_that("Character limitation works correctly", {
  distance_scores <- Matrix::Matrix(runif(1000, min=0, max=1), 100, 100)
  char_lim <- 10
  p <- distance_heatmap(distance_scores, chars_limit = char_lim)
  expect_type(p, "list")
})
