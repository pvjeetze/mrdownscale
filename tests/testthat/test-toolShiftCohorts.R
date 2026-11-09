test_that("cohorts move on and the last cohort accumulates", {
  stock <- cohortHr(1, 2, 3)
  expect_equal(toolShiftCohorts(stock, 0), stock)
  expect_equal(as.vector(toolShiftCohorts(stock, 1)), c(0, 1, 5))
  expect_equal(as.vector(toolShiftCohorts(stock, 2)), c(0, 0, 6))
  expect_equal(as.vector(toolShiftCohorts(stock, 5)), c(0, 0, 6))
  expect_identical(dimnames(toolShiftCohorts(stock, 1)), dimnames(stock))
})
