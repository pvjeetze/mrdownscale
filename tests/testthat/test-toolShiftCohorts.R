test_that("cohorts move on and the last cohort accumulates", {
  stock <- cohortStock(1, 2, 3)
  expect_equal(toolShiftCohorts(stock, 0), stock)
  expect_equal(unname(toolShiftCohorts(stock, 1)[1, ]), c(0, 1, 5))
  expect_equal(unname(toolShiftCohorts(stock, 2)[1, ]), c(0, 0, 6))
  expect_equal(unname(toolShiftCohorts(stock, 5)[1, ]), c(0, 0, 6))
})
