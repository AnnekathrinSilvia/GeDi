data(macrophage_topGO_example_small, package = "GeDi")
data(ppi_macrophage_topGO_example_small, package = "GeDi")
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

# testServer(gedi_server, {
#   expect_true(is.null(input$inputgenesetfile))
#   expect_true(is.null(input$bin_gs_hist[[1]]))
#   expect_true(is.null(input$bin_gs_hist[[2]]))
#   expect_true(is.null(input$bindwidth_hist))
#   expect_true(is.null(input$species))
#   expect_true(is.null(input$scoringmethod))
#   expect_true(is.null(input$cluster_method_dendro))
#   expect_true(is.null(input$similarityScores))
#   expect_true(is.null(input$select_clustering))
#   expect_true(is.null(input$graphColoring))
#   expect_true(is.null(input$alt_names_start))
#   expect_true(is.null(input$alt_name_genesets))
#   expect_true(is.null(input$alt_name_genes))
#   expect_true(is.null(input$filter_genesets))
#   expect_true(is.null(input$select_filter_genesets_threshold))
# })
#
# test_server <- testServer(gedi_server, function(input){
#     session$setInputs(btn_info_session = 1)
#     expect_true(input$btn_info_session == 1)
# })

