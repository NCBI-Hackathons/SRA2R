library(SRA2R)
context("getPileUp")

test_that("getPileUp returns pileups from a run", {
  expect_equal(class(getPileUp('SRR1607152', '20', 1000000, 1000100)),"data.frame")
  expect_equal(nrow(getPileUp('SRR1607152', '20', 1000000, 1000100)),101)
  expect_equal(ncol(getPileUp('SRR1607152', '20', 1000000, 1000100)),6)
})
