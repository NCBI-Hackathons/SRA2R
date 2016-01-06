library(SRA2R)
context("getReference")

test_that("getReference returns references from a run", {
  expect_equal(class(getReference('SRR390728')),"data.frame")
  expect_equal(nrow(getReference('SRR390728')),25)
  expect_equal(ncol(getReference('SRR390728')),4)
  expect_equal(nrow(getReference('SRR1607152')),86)
  expect_equal(ncol(getReference('SRR1607152')),4)
})
