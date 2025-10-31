test_that("distanceDendro-no scores", {
  expect_error(distanceDendro(NULL))
  expect_error(distanceDendro(list()))
})

test_that("distanceDendro runs correctly", {
  data(scores_macrophage_topGO_example_small, package = "GeDi")
  dendro <- distanceDendro(scores_macrophage_topGO_example_small,
    cluster_method = "single"
  )
  expect_type(dendro, "list")
  dendro <- distanceDendro(scores_macrophage_topGO_example_small,
    cluster_method = "average"
  )
  expect_type(dendro, "list")
  expect_error(distanceDendro(scores_macrophage_topGO_example_small,
    cluster_method = "test"
  ))
})
