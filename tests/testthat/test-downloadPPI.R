test_that("no id", {
  expect_true(is.na(getId("")))
})

test_that("Homo sapiens ID", {
  expect_true(getId("Homo sapiens") == "9606")
})

test_that("Homo Sapiens Misspelled", {
  expect_true(is.na(getId("Homo Spasiens")))
})


test_that("StringDB of no species", {
  expect_error(getStringDB(NA))
  expect_error(getStringDB("NA"))
})


test_that("PPI retrieval works", {
  stringdb <- getStringDB(9606)
  stringdb

  anno_df <- getAnnotation(stringdb)

  expect_type(anno_df, "list")
  genes <- c(c("CFTR", "RALA"), c("CACNG3", "ITGA3"), c("DVL2"))
  ppi <- getPPI(genes, stringdb, anno_df)

  expect_type(ppi, "list")
})
