fillingCase <- function() {
  # unit A with three grid cells, unit B with one; two items
  mapping <- data.frame(cell = c("c1", "c2", "c3", "c4"), cluster = c("A", "A", "A", "B"))
  prior <- new.magpie(mapping[["cell"]], NULL, c("x", "y"), fill = c(0.2, 0.1, 0, 1, 1, 1, 1, 0))
  cap <- new.magpie(mapping[["cell"]], NULL, c("x", "y"), fill = 1)
  return(list(mapping = mapping, prior = prior, cap = cap))
}

test_that("totals are distributed in proportion to the prior, within the caps", {
  x <- fillingCase()
  total <- new.magpie(c("A", "B"), NULL, c("x", "y"), fill = c(0.6, 0.5, 1.5, 0))
  out <- toolAllocateByFilling(total, x[["prior"]], x[["cap"]], x[["mapping"]])
  # A, x: 0.6 in proportion 2:1
  expect_equal(as.vector(out[c("c1", "c2", "c3"), , "x"]), c(0.4, 0.2, 0))
  # A, y: equal priors, 0.5 each
  expect_equal(as.vector(out[c("c1", "c2", "c3"), , "y"]), c(0.5, 0.5, 0.5))
  # B: prior 1 for x; zero total for y
  expect_equal(as.vector(out["c4", , ]), c(0.5, 0))
  expect_identical(dimnames(out), dimnames(x[["prior"]]))
})

test_that("grid cells with a prior fill up first, the rest goes to the others by their remaining cap", {
  x <- fillingCase()
  total <- new.magpie(c("A", "B"), NULL, "x", fill = c(2.5, 0))
  out <- toolAllocateByFilling(total, x[["prior"]][, , "x"], x[["cap"]][, , "x"], x[["mapping"]])
  # c1 and c2 reach their cap of 1, the rest of 0.5 goes to c3, which has no prior
  expect_equal(as.vector(out[c("c1", "c2", "c3"), , ]), c(1, 1, 0.5))
  # beyond all caps: every grid cell at its cap
  total["A", , ] <- 5
  out <- toolAllocateByFilling(total, x[["prior"]][, , "x"], x[["cap"]][, , "x"], x[["mapping"]])
  expect_equal(as.vector(out[c("c1", "c2", "c3"), , ]), c(1, 1, 1))
})

test_that("years are filled independently", {
  x <- fillingCase()
  prior <- mbind(setYears(x[["prior"]][, , "x"], "y2000"), setYears(x[["prior"]][, , "x"], "y2005"))
  cap <- prior * 0 + 1
  total <- new.magpie(c("A", "B"), c("y2000", "y2005"), "x", fill = c(0.6, 0.5, 2.5, 0))
  out <- toolAllocateByFilling(total, prior, cap, x[["mapping"]])
  expect_equal(as.vector(out[c("c1", "c2", "c3"), "y2000", ]), c(0.4, 0.2, 0))
  expect_equal(as.vector(out[c("c1", "c2", "c3"), "y2005", ]), c(1, 1, 0.5))
})

test_that("the result equals toolAllocateByCrossEntropy with one hard cap per grid cell", {
  x <- fillingCase()
  total <- new.magpie(c("A", "B"), NULL, "x", fill = c(1.1, 0.5))
  prior <- x[["prior"]][, , "x"]
  cap <- x[["cap"]][, , "x"]
  cap["c1", , ] <- 0.5
  closed <- toolAllocateByFilling(total, prior, cap, x[["mapping"]])
  iterative <- toolAllocateByCrossEntropy(total, prior, x[["mapping"]], caps = list(list(group = "x", cap = cap)),
                                          tol = 1e-12)[["flowsHr"]]
  expect_equal(closed, iterative, tolerance = 1e-9)
})
