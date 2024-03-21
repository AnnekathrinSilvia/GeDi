data(macrophage_topGO_example_small, package = "GeDi")
data(ppi_macrophage_topGO_example_small, package = "GeDi")
data(scores_macrophage_topGO_example_small, package = "GeDi")
library("GeDi")

test_that("Shiny app is generated", {
  expect_s3_class(GeDi(), "shiny.appobj")
})

test_that("Shiny app is generated with input", {
  expect_s3_class(GeDi(genesets = macrophage_topGO_example_small), "shiny.appobj")

  expect_s3_class(GeDi(genesets = macrophage_topGO_example_small,
                       ppi = ppi_macrophage_topGO_example_small), "shiny.appobj")

  expect_s3_class(GeDi(genesets = macrophage_topGO_example_small,
                       ppi = ppi_macrophage_topGO_example_small,
                       distance_scores = scores_macrophage_topGO_example_small), "shiny.appobj")
})

