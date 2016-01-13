library(SRA2R)
context("getFastqReads")

test_that("getFastqReads returns reads from a run", {
  #expect the need to update when object class is updated
  expect_equal(class(getFastqReads('SRR1607152',10)),"list")
  expect_equal(listLen(getFastqReads('SRR1607152',10)),20)
})
