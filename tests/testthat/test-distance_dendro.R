test_that("distance_dendro-no scores", {
  expect_error(distance_dendro(NULL))
  expect_error(distance_dendro(list()))
})

test_that("distance_dendro runs correctly", {
  data(scores_macrophage_topGO_example_small, package = "GeDi")
  dendro <- distance_dendro(scores_macrophage_topGO_example_small,
    cluster_method = "single"
  )
  expect_type(dendro, "list")
  dendro <- distance_dendro(scores_macrophage_topGO_example_small,
    cluster_method = "average"
  )
  expect_type(dendro, "list")
  expect_error(distance_dendro(scores_macrophage_topGO_example_small,
    cluster_method = "test"
  ))
})
