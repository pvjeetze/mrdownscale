test_that("along the transitions without harmonisation, the land is summed over the sub-types", {
  landHr <- new.magpie(c("1.1.AAA", "2.2.AAA"), c("y1985", "y2000", "y2005"),
                       c("crop.total", "past.total", "forestry.aff", "forestry.plant", "primforest.total",
                         "secdforest.total", "urban.total", "other.othernat", "other.youngsecdf"), fill = 0.1)
  localMagpieSources(list(), outputs = list("LandTransitionsHighRes.input=magpie" = list(landHr = landHr)))
  out <- calcLandHighRes("magpie", "landuseinit", NULL, 2005, "none", "transitions")[["x"]]
  # the start (y1985) stays with the selected years
  expect_equal(getYears(out), c("y1985", "y2005"))
  expect_equal(getNames(out), c("crop", "past", "forestry", "primforest", "secdforest", "urban", "other"))
  expect_equal(as.vector(out[, , "forestry"]), rep(0.2, 4))
  expect_error(calcLandHighRes("magpie", "luh3", NULL, NULL, "none", "transitions"), "supports only target")
})
