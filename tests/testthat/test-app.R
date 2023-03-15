data(macrophage_topGO_example_small, package = "GeDi")
data(ppi_macrophage_topGO_example_small, package = "GeDi")

test_that("Shiny app is generated", {
  expect_type(
    GeDi(),
    "list"
  )
})

test_that("Shiny app is generated with input",{
  expect_type(GeDi(genesets = macrophage_topGO_example_small),
              "list")
  expect_type(GeDi(genesets = macrophage_topGO_example_small,
                   ppi = ppi),
              "list")
})






