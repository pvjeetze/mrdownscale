test_that("age classes follow the grid cell land and the age classes the model removes", {
  case <- ageCase()
  localMagpieSources(case[["sources"]], outputs = list("LandTransitionsHighRes.input=magpie" = case[["land"]]))
  out <- calcAgeClassesHighRes("magpie")[["x"]]
  expect_equal(names(out), c("secdforest", "forestry", "other", "diagnostics"))
  # grid cell 1 loses the oldest age class (acx after the shift), grid cell 2 keeps all
  expect_equal(as.vector(out[["secdforest"]]["c1", "y2005", ]), c(0, 0, 0.5, 0), tolerance = 1e-8)
  expect_equal(as.vector(out[["secdforest"]]["c2", "y2005", ]), c(0, 0, 0.5, 0.5), tolerance = 1e-8)
  # new other land enters the youngest age class
  expect_equal(as.vector(out[["other"]]["c2", "y2005", ]), c(0.2, 0, 0, 0), tolerance = 1e-8)
  diagnostics <- out[["diagnostics"]]
  expect_equal(unique(diagnostics[["type"]]), c("secdforest", "plant", "othernat", "youngsecdf"))
  expect_lt(max(diagnostics[["gridCellDeviation"]]), 1e-10)
  expect_lt(max(diagnostics[["clusterDeviation"]]), 1e-8)
})
