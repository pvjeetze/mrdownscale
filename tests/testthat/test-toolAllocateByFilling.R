fillingCase <- function() {
  # unit A with three grid cells, unit B with one; two columns
  mapping <- data.frame(cell = c("c1", "c2", "c3", "c4"), cluster = c("A", "A", "A", "B"))
  prior <- matrix(c(0.2, 0.1, 0, 1, 1, 1, 1, 0), 4, 2, dimnames = list(mapping[["cell"]], c("x", "y")))
  cap <- matrix(1, 4, 2, dimnames = dimnames(prior))
  return(list(mapping = mapping, prior = prior, cap = cap))
}

test_that("totals are distributed in proportion to the prior, within the caps", {
  x <- fillingCase()
  total <- matrix(c(0.6, 0.5, 1.5, 0), 2, 2, dimnames = list(c("A", "B"), c("x", "y")))
  out <- toolAllocateByFilling(total, x[["prior"]], x[["cap"]], x[["mapping"]])
  # A, x: 0.6 in proportion 2:1
  expect_equal(unname(out[1:3, "x"]), c(0.4, 0.2, 0))
  # A, y: equal priors, 0.5 each
  expect_equal(unname(out[1:3, "y"]), c(0.5, 0.5, 0.5))
  # B: prior 1 for x; zero total for y
  expect_equal(unname(out[4, ]), c(0.5, 0))
  expect_equal(dimnames(out), dimnames(x[["prior"]]))
})

test_that("grid cells with a prior fill up first, the rest goes to the others by their remaining cap", {
  x <- fillingCase()
  total <- matrix(c(2.5, 0), 2, 1, dimnames = list(c("A", "B"), "x"))
  out <- toolAllocateByFilling(total, x[["prior"]][, "x", drop = FALSE], x[["cap"]][, "x", drop = FALSE],
                               x[["mapping"]])
  # c1 and c2 reach their cap of 1, the rest of 0.5 goes to c3, which has no prior
  expect_equal(unname(out[1:3, "x"]), c(1, 1, 0.5))
  # beyond all caps: every grid cell at its cap
  total[1, 1] <- 5
  out <- toolAllocateByFilling(total, x[["prior"]][, "x", drop = FALSE], x[["cap"]][, "x", drop = FALSE],
                               x[["mapping"]])
  expect_equal(unname(out[1:3, "x"]), c(1, 1, 1))
})

test_that("the result equals toolAllocateByCrossEntropy with one hard cap per grid cell", {
  x <- fillingCase()
  total <- matrix(c(1.1, 0.5), 2, 1, dimnames = list(c("A", "B"), "x"))
  prior <- x[["prior"]][, "x", drop = FALSE]
  cap <- x[["cap"]][, "x", drop = FALSE]
  cap["c1", "x"] <- 0.5
  closed <- toolAllocateByFilling(total, prior, cap, x[["mapping"]])
  iterative <- toolAllocateByCrossEntropy(total, prior, x[["mapping"]], caps = list(list(group = "x", cap = cap)),
                                          tol = 1e-12)[["flowsHr"]]
  expect_equal(closed, iterative[rownames(closed), , drop = FALSE], tolerance = 1e-9)
})
