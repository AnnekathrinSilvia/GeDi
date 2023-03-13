test_that("No genesets - getGenes", {
  genesets <- list()
  expect_true(is.null(getGenes(genesets = genesets)))
})

test_that("No column named Genes - getGenes", {
  df <- data.frame(
    Geneset = c(
      "Cell Cycle",
      "Biological Process",
      "Mitosis"
    ),
    Gen = c(
      c("PDHB,VARS2,IARS2"),
      c("LARS,LARS2"),
      c("IARS,SUV3")
    )
  )
  expect_error(getGenes(genesets = df))
  genes <- getGenes(genesets = df, gene_name = "Gen")
  expect_equal(length(genes), 3)
  expect_type(genes, "list")
})

test_that("getGenes runs correctly", {
  df <- data.frame(
    Geneset = c(
      "Cell Cycle",
      "Biological Process",
      "Mitosis"
    ),
    Genes = c(
      c("PDHB,VARS2,IARS2"),
      c("LARS,LARS2"),
      c("IARS,SUV3")
    )
  )
  genes <- getGenes(genesets = df)
  expect_equal(length(genes), 3)
  expect_type(genes, "list")
})

test_that("Remove no genesets", {
  df <- data.frame(
    Geneset = c(
      "Cell Cycle",
      "Biological Process",
      "Mitosis"
    ),
    Genes = c(
      c("PDHB,VARS2,IARS2"),
      c("LARS,LARS2"),
      c("IARS,SUV3")
    )
  )
  test_df <- .filterGenesets(remove = list(), df)
  expect_identical(df, test_df[[1]])
  test_df <- .filterGenesets(remove = "ABC CDF", df)
  expect_identical(df, test_df[[1]])
})

test_that(".filterGenesets runs correctly", {
  df <- data.frame(
    Geneset = c(
      "Cell Cycle",
      "Biological Process",
      "Mitosis"
    ),
    Genes = c(
      c("PDHB,VARS2,IARS2"),
      c("LARS,LARS2"),
      c("IARS,SUV3")
    )
  )
  genes <- getGenes(df)
  genes <- genes[c(-3)]
  test_df <- .filterGenesets(remove = "Mitosis", df)
  expect_identical(c("Cell Cycle", "Biological Process"), test_df[[2]])
  expect_identical(genes, test_df[[3]])
})

test_that(".getNumberCores does not take all cores", {
  available_cores <- parallel::detectCores()
  overhead <- available_cores + 100

  cores <- .getNumberCores(n_cores = overhead)
  expect_lte(cores, available_cores)
  cores <- .getNumberCores()
  expect_lte(cores, available_cores)
})

test_that(".getNumberCores returns correct number", {
  expect_equal(.getNumberCores(n_cores = 3), 3)
})

test_that(".sepGuesser returns correct separator", {
  sep <- .sepguesser(system.file("extdata", "intro_data_upload.txt", package = "GeDi"))
  expect_equal(sep, ";")
})

test_that("checkPPI", {
  expect_error(.checkPPI(list()))
  expect_error(.checkPPI())
  df <- data.frame(c(1, 2, 3), c("a", "b", "c"))
  expect_error(.checkPPI(df))
  df <- data.frame(c(1, 2, 3), c("a", "b", "c"), c("d", "e", "f"))
  expect_error(.checkPPI(df))
  df <- data.frame(Gene1 = c("a", "b", "c"),
                   Gene2 = c("a", "b", "c"),
                   combined_score = c(1, 2, 3))
  expect_identical(.checkPPI(df), df)
})

test_that("checkGenesets", {
  expect_error(.checkGenesets(list()))
  expect_error(.checkGenesets())
  df <- data.frame(c("a", "b"), c("c", "d"))
  expect_error(.checkGenesets(df))
  df <- data.frame(Geneset = c(1, 2), Genes = c("a", "b"))
  expect_error(.checkGenesets(df))
  df <- data.frame(Genes = c(1, 2), Geneset = c("a", "b"))
  expect_error(.checkGenesets(df))
  df <- data.frame(
    Geneset = c(
      "Cell Cycle",
      "Biological Process",
      "Mitosis"
    ),
    Genes = c(
      c("PDHB,VARS2,IARS2"),
      c("LARS,LARS2"),
      c("IARS,SUV3")
    )
  )
  expect_identical(df, .checkGenesets(df))
})
