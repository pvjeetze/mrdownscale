test_that("cropland is bounded by the available cropland reduced by the SNV share and its fader", {
  case <- magpieCase()
  sources <- magpieSources(case, availableCropland = c(1, 2, 3))
  sources[["snv"]][["shareIso"]][["A"]] <- 0.2
  sources[["snv"]][["fader"]][, "y2005", ] <- 0.5
  localMagpieSources(sources)
  constraints <- calcLandConstraints("magpie")[["x"]]
  expect_equal(names(constraints), c("potentialForest", "crop"))
  expect_equal(as.vector(constraints[["crop"]][["upper"]][, "y2000", ]), c(1, 2, 3))
  expect_equal(as.vector(constraints[["crop"]][["upper"]][, "y2005", ]), c(1, 2, 3) * 0.9)
  expect_equal(constraints[["crop"]][["priority"]], 3)
})

test_that("potential forest without climate change keeps the values of 1995", {
  case <- magpieCase()
  sources <- magpieSources(case, config = list(c35_pot_forest_scenario = "nocc"))
  cells <- case[["mapping"]][["cell"]]
  sources[["potentialForest"]] <- new.magpie(cells, c("y1995", "y2000", "y2005"), NULL, fill = rep(1:3, each = 3))
  localMagpieSources(sources)
  upper <- calcLandConstraints("magpie")[["x"]][["potentialForest"]][["upper"]]
  expect_equal(getYears(upper), c("y1995", "y2000", "y2005"))
  expect_equal(as.vector(upper), rep(1, 9))
})

test_that("protected natural vegetation and pasture are lower bounds where the run has conservation data", {
  case <- magpieCase()
  sources <- magpieSources(case)
  sources[["available"]] <- c(sources[["available"]], "wdpa")
  cells <- case[["mapping"]][["cell"]]
  types <- c("crop", "past", "forestry", "primforest", "secdforest", "urban", "other")
  protection <- new.magpie(cells, c("y2000", "y2005"), types, fill = 0)
  protection[, , "primforest"] <- 0.1
  protection[, , "secdforest"] <- 0.2
  protection[, , "other"] <- 0.3
  protection[, , "past"] <- 0.4
  localMagpieSources(sources, outputs = list("LandConservationHighRes.input=magpie" = protection))
  constraints <- calcLandConstraints("magpie")[["x"]]
  expect_equal(constraints[["natveg"]][["types"]], c("primforest", "secdforest", "other"))
  expect_equal(as.vector(constraints[["natveg"]][["lower"]]), rep(0.6, 6))
  expect_equal(as.vector(constraints[["past"]][["lower"]]), rep(0.4, 6))
  # protection is exceeded first, then potential forest, then cropland
  expect_equal(vapply(constraints, function(x) x[["priority"]], numeric(1)),
               c(potentialForest = 2, crop = 3, natveg = 1, past = 1))
})
