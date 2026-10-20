test_that("without cell totals, cohort flows follow the cell stock of the cohort", {
  stock <- cohortStock(1, 0, 2,
                       3, 1, 2,
                       1, 1, 1)
  out <- toolAllocateCohortFlows(stock, cohortFlow(2, 0.5, 1, 0.5, 0.5, 0.5), cellMapping(c("A", "A", "B")))
  expect_equal(unname(out[1:2, "ac0"]), c(0.5, 1.5))
  expect_equal(unname(out[1:2, "ac10"]), c(0.5, 0.5))
  expect_equal(unname(colSums(out[1:2, ])), c(2, 0.5, 1))
  expect_equal(unname(out[3, ]), c(0.5, 0.5, 0.5))
  # a cluster flow above the cluster stock of the cohort is not allocated
  out <- toolAllocateCohortFlows(stock, cohortFlow(5, 0, 0, 0, 0, 0), cellMapping(c("A", "A", "B")))
  expect_equal(unname(out[1:2, "ac0"]), c(1, 3))
})

test_that("with cell totals, both margins hold where the stocks allow it", {
  stock <- cohortStock(1, 1, 0,
                       1, 1, 2,
                       0, 2, 2)
  flowLr <- cohortFlow(1, 1, 2)
  totalHr <- c(1, 1.5, 1.5)
  out <- toolAllocateCohortFlows(stock, flowLr, cellMapping(rep("A", 3)), totalHr)
  expect_equal(unname(rowSums(out)), totalHr, tolerance = 1e-8)
  expect_equal(colSums(out), flowLr[1, ], tolerance = 1e-6)
  expect_true(all(out <= stock + 1e-10))
})

test_that("cell totals hold where the stocks do not allow the cluster flows of the cohorts", {
  stock <- cohortStock(1, 0, 0,
                       0, 0, 1)
  # the cluster removes only ac10, but cell 1 has none: its total comes from ac0
  out <- toolAllocateCohortFlows(stock, cohortFlow(0, 0, 1), cellMapping(c("A", "A")), totalHr = c(0.5, 0.5))
  expect_equal(unname(rowSums(out)), c(0.5, 0.5), tolerance = 1e-8)
  expect_equal(unname(out[1, ]), c(0.5, 0, 0), tolerance = 1e-8)
  expect_true(all(out <= stock + 1e-10))
})

test_that("cell totals without cluster flows follow the cohort composition of the stock", {
  stock <- cohortStock(1, 3, 0,
                       1, 1, 0)
  out <- toolAllocateCohortFlows(stock, cohortFlow(0, 0, 0), cellMapping(c("A", "A")), totalHr = c(1, 0))
  expect_equal(unname(rowSums(out)), c(1, 0), tolerance = 1e-8)
  expect_equal(unname(out[1, ]), c(0.25, 0.75, 0), tolerance = 1e-6)
})
