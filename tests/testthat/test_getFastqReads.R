library(SRA2R)
context("getFastqReads")

test_that("getFastqReads returns reads from a run", {
  expect_equal(getFastqReads('SRR1607152',10),20)
})
