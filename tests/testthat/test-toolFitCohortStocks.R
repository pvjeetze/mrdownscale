test_that("cell stocks are fitted to the coarse stocks, keeping the cell totals", {
  stock <- cohortHr(1, 0, 0,
                    0, 1, 1,
                    2, 0, 0)
  stockLr <- cohortLr(1, 1, 1,
                      0, 1, 1)
  out <- toolFitCohortStocks(stock, stockLr, cellMapping(c("A", "A", "B")), tol = 1e-10)
  expect_equal(as.vector(dimSums(out, dim = 3)), c(1, 2, 2), tolerance = 1e-12)
  expect_equal(as.vector(dimSums(out[c("c1", "c2"), , ], dim = 1)), as.vector(stockLr["A", , ]), tolerance = 1e-8)
  # cell 3 holds only ac0, but its unit has only ac5 and ac10: the seed lets it take them
  expect_equal(as.vector(out["c3", , ]), c(0, 1, 1), tolerance = 1e-6)
})

test_that("rounding residues in the cell stocks give the result of exact zeros", {
  # cell 1 holds only ac0, which its unit lacks: the seed lets it take ac5 and ac10 as cell 2
  stock <- cohortHr(2, 0, 0,
                    0, 1, 1)
  stockLr <- cohortLr(0, 2, 2)
  mapping <- cellMapping(c("A", "A"))
  out <- toolFitCohortStocks(stock, stockLr, mapping, tol = 1e-10)
  expect_equal(as.vector(out), c(0, 0, 1, 1, 1, 1), tolerance = 1e-6)
  # a residue of ac5 in cell 1 must not keep it from the seed (it would take all of its total as ac5)
  residue <- stock
  residue["c1", , "ac5"] <- 1e-17
  expect_equal(as.vector(toolFitCohortStocks(residue, stockLr, mapping, tol = 1e-10)), as.vector(out),
               tolerance = 1e-12)
})
