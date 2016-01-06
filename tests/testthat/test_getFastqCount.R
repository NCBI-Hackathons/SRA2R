library(SRA2R)
context("getFastqCount")

test_that("getFastqCount returns a count for the number of reads in a run", {
  expect_equal(getFastqCount('SRR000123'),4583)
  expect_equal(getFastqCount('SRR1607152'),78377869)
})
