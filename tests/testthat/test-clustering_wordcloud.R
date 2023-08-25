data(macrophage_topGO_example,
     package = "GeDi",
     envir = environment())

test_that("No genesets are provided", {
  expect_error(enrichmentWordcloud(NULL))
})

test_that("EnrichmentWordcloud is correctly created", {
  p <- enrichmentWordcloud(macrophage_topGO_example)
  expect_type(p, "list")

  data <- macrophage_topGO_example
  colnames(data) <- c("Genesets", "Description", "Annotated", "Significant", "Expected",
                      "Rank", "p_value", "padj", "Genes")
  p <- enrichmentWordcloud(data)
  expect_type(p, "list")
})

test_that("EnrichmentWordcloud works correctly if only Genesets and Genes column are given", {
  data <- macrophage_topGO_example[, c(1, 9)]
  rownames(data) <- data$Genesets
  p <- enrichmentWordcloud(data)
  expect_type(p, "list")
})
