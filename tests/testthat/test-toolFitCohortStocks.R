test_that("cell stocks are fitted to the coarse stocks, keeping the cell totals", {
  stock <- cohortStock(1, 0, 0,
                       0, 1, 1,
                       2, 0, 0)
  stockLr <- cohortFlow(1, 1, 1,
                        0, 1, 1)
  out <- toolFitCohortStocks(stock, stockLr, cellMapping(c("A", "A", "B")), tol = 1e-10)
  expect_equal(as.vector(dimSums(out, dim = 3)), c(1, 2, 2), tolerance = 1e-12)
  expect_equal(as.vector(dimSums(out[c("c1", "c2"), , ], dim = 1)), as.vector(stockLr["A", , ]), tolerance = 1e-8)
  # cell 3 holds only ac0, but its unit has only ac5 and ac10: the seed lets it take them
  expect_equal(as.vector(out["c3", , ]), c(0, 1, 1), tolerance = 1e-6)
})
