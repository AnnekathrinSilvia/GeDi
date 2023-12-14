data(scores_macrophage_topGO_example_small, package = "GeDi")
test_that("No distance scores", {
  expect_error(distanceHeatmap(NULL))
})

test_that("Heatmap is correctly created", {
  p <- distanceHeatmap(scores_macrophage_topGO_example_small)
  expect_type(p, "S4")
})

test_that("Character limitation works correctly", {
  p <- distanceHeatmap(scores_macrophage_topGO_example_small,
    chars_limit = 5
  )
  expect_type(p, "S4")
})
