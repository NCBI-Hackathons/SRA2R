library(SRA2R)
context("getReference")

test_that("getReference returns references from a run", {
  expect_equal(nrow(getReference('SRR390728')),25)
})
