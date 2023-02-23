test_that("No distance scores", {
  distance_scores <- NULL
  expect_error(distance_heatmap(distance_scores = distance_scores))
})

test_that("Heatmap is correctly created", {
  names <- c(1:10)
  names <- as.character(names)
  distance_scores <- Matrix::Matrix(stats::runif(100, min = 0, max = 1), 10, 10)
  rownames(distance_scores) <- colnames(distance_scores) <- names
  p <- distance_heatmap(distance_scores)
  expect_type(p, "list")
})

test_that("Character limitation works correctly", {
  names <- c(1:10)
  names <- as.character(names)
  distance_scores <- Matrix::Matrix(stats::runif(100, min = 0, max = 1), 10, 10)
  rownames(distance_scores) <- colnames(distance_scores) <- names
  char_lim <- 10
  p <- distance_heatmap(distance_scores, chars_limit = char_lim)
  expect_type(p, "list")
})
