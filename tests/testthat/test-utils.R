test_that("No genesets - getGenes", {
  genesets <- list()
  expect_equal(length(getGenes(genesets = genesets)), 0)
})

test_that("No column named Genes - getGenes", {
  df <- data.frame(
    Geneset = c("Cell Cycle",
                "Biological Process",
                "Mitosis"),
    Gen =  c(c("PDHB,VARS2,IARS2"),
             c("LARS,LARS2"),
             c("IARS,SUV3"))
  )
  expect_error(getGenes(genesets = df))
})

test_that("getGenes runs correctly", {
  df <- data.frame(
    Geneset = c("Cell Cycle",
                "Biological Process",
                "Mitosis"),
    Genes =  c(c("PDHB,VARS2,IARS2"),
               c("LARS,LARS2"),
               c("IARS,SUV3"))
  )
  genes <- getGenes(genesets = df)
  expect_equal(length(genes), 3)
  expect_type(genes, "list")
})




