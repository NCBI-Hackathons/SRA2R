library(SRA2R)
context("getReads")

test_that("getReads returns reads from a run", {
  expect_equal(getFastqCount('SRR000123'),4583)
})
