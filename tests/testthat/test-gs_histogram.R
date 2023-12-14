test_that("gsHistogram-no genes", {
  genes <- NULL
  expect_error(gsHistogram(genes, gs_names = list()))
  expect_error(gsHistogram(list(), list()))
})

test_that("gsHistogram-no geneset names", {
  genes <- list(
    c("PDHB", "VARS2"), c("IARS2", "PDHA1"),
    c("AAAS", "ABCE1"), c("ABI1", "AAR2"), c("AATF", "AMFR"),
    c("BMS1", "DAP3"), c("AURKAIP1", "CHCHD1"), c("IARS2"),
    c("AHI1", "ALMS1")
  )
  expect_error(gsHistogram(genes, gs_names = list()))
  expect_error(gsHistogram(genes, gs_names = c("a", "b")))
})

test_that("gsHistogram runs correctly", {
  gs_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
  genes <- list(
    c("PDHB", "VARS2"), c("IARS2", "PDHA1"),
    c("AAAS", "ABCE1"), c("ABI1", "AAR2"), c("AATF", "AMFR"),
    c("BMS1", "DAP3"), c("AURKAIP1", "CHCHD1"), c("IARS2"),
    c("AHI1", "ALMS1")
  )
  histogram <- gsHistogram(genes, gs_names = gs_names)
  expect_type(histogram, "list")

  histogram <- gsHistogram(genes, gs_names, start = 1, end = 5)
  expect_type(histogram, "list")
})
