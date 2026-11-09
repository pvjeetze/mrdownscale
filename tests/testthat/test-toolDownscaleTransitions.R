test_that("downscaling along the net transitions keeps the unit states and the grid cell area", {
  # two units of two grid cells; unit A loses 0.4 crop to secdf from 2020 to 2025, unit B gains 0.2 crop
  types <- c("crop", "secdf", "urban")
  x <- new.magpie(c("A.1", "B.2"), c("y2020", "y2025", "y2030"), types, fill = 0)
  x["A.1", , ] <- c(1, 0.6, 0.6, 1, 1.4, 1.4, 0, 0, 0)
  x["B.2", , ] <- c(0.5, 0.7, 0.7, 1.5, 1.3, 1.3, 0, 0, 0)
  cells <- c("1p25.1p25", "2p25.2p25", "3p25.3p25", "4p25.4p25")
  xTarget <- new.magpie(cells, "y2020", types, fill = 0)
  xTarget[, , "crop"] <- c(0.8, 0.2, 0.5, 0)
  xTarget[, , "secdf"] <- c(0.2, 0.8, 0.5, 1)
  mapping <- data.frame(cell = cells, lowRes = c("A.1", "A.1", "B.2", "B.2"))

  out <- toolDownscaleTransitions(x, xTarget, mapping)
  expect_equal(getYears(out), c("y2025", "y2030"))
  expect_equal(unname(getSets(out)), c("x", "y", "year", "data"))
  # the unit states hold, and the area of each grid cell stays the same
  statesLr <- toolAggregate(out, data.frame(from = cells, to = mapping[["lowRes"]]), from = "from", to = "to")
  expect_equal(as.vector(statesLr[c("A.1", "B.2"), , ]), as.vector(x[, -1, ]), tolerance = 1e-6)
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(1, 8), tolerance = 1e-8)
  # unit A loses crop in proportion to the crop of its grid cells; B gains it from secdf
  expect_equal(as.vector(out[c("1p25.1p25", "2p25.2p25"), "y2025", "crop"]), c(0.48, 0.12), tolerance = 1e-6)
  expect_true(all(out >= 0))
})
