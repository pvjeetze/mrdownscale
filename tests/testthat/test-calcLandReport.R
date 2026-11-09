# Land of calcLandTransitionsHighRes for two grid cells: cropland 1 and 2, crop area 60 %, fallow 30 %,
# tree cover 10 %; crop area 2/3 rainfed food crops, 1/6 irrigated food crops, 1/6 rainfed bioenergy crops.
transitionsCase <- function() {
  cells <- c("1.1.AAA", "2.2.AAA")
  years <- c("y1985", "y2000")
  landHr <- new.magpie(cells, years, c("crop.total", "past.total", "forestry.aff", "primforest.total",
                                       "secdforest.total", "urban.total", "other.othernat"), fill = 1)
  landHr[, , "crop.total"] <- c(1, 2, 1, 2)
  cropSplitHr <- new.magpie(cells, years, c("crop_area", "crop_fallow", "crop_treecover"), fill = 0)
  for (i in 1:3) {
    cropSplitHr[, , i] <- dimSums(landHr[, , "crop"], dim = 3) * c(0.6, 0.3, 0.1)[i]
  }
  cropAreaHr <- new.magpie(cells, years, c("tece.rainfed", "tece.irrigated", "betr.rainfed", "betr.irrigated",
                                           "begr.rainfed", "begr.irrigated"), fill = 0)
  shares <- c(2 / 3, 1 / 6, 1 / 6, 0, 0, 0)
  for (i in 1:6) {
    cropAreaHr[, , i] <- cropSplitHr[, , "crop_area"] * shares[i]
  }
  return(list(landHr = landHr, splitHr = list(crop = cropSplitHr, crop_area = cropAreaHr)))
}

test_that("the land split groups the crop area into food and bioenergy crops and covers the land", {
  localMagpieSources(list(), outputs = list("LandTransitionsHighRes.input=magpie" = transitionsCase()))
  split <- calcLandReport("downscaledmagpieSplit", "magpie", NULL, NULL, "none", "transitions")[["x"]]
  expect_equal(getNames(split),
               c("crop.kfo.rainfed", "crop.kfo.irrigated", "crop.kbe.rainfed", "crop.kbe.irrigated",
                 "crop.fallow.total", "crop.treecover.total", "past.total.total", "forestry.aff.total",
                 "primforest.total.total", "secdforest.total.total", "urban.total.total", "other.othernat.total"))
  expect_equal(getSets(split, fulldim = FALSE)[3], "land.type.water")
  # grid cell 2: cropland 2
  expect_equal(as.vector(split["2.2.AAA", "y2000", ]), c(0.8, 0.2, 0.2, 0, 0.6, 0.2, rep(1, 6)))
})

test_that("crop area is given as share of the grid cell area, in the years of the model", {
  localMagpieSources(list(), outputs = list("LandTransitionsHighRes.input=magpie" = transitionsCase()))
  share <- calcLandReport("downscaledmagpieCroparea", "magpie", NULL, NULL, "none", "transitions")[["x"]]
  expect_equal(getYears(share), "y2000")
  expect_equal(getSets(share, fulldim = FALSE)[3], "kcr.w")
  expect_equal(as.vector(share["2.2.AAA", , "tece.rainfed"]), 0.8 / 8)
})

test_that("the land pools of the run keep its start and grid cells", {
  localMagpieSources(list(), outputs = list("LandTransitionsHighRes.input=magpie" = transitionsCase()))
  pools <- calcLandReport("downscaledmagpie", "magpie", c(2015, 2050), NULL, "none", "transitions")[["x"]]
  expect_equal(getYears(pools), c("y1985", "y2000"))
  expect_equal(getCells(pools), c("1.1.AAA", "2.2.AAA"))
  expect_equal(getNames(pools), c("crop", "past", "forestry", "primforest", "secdforest", "urban", "other"))
  expect_equal(as.vector(pools[, , "crop"]), c(1, 2, 1, 2))
})

test_that("the MAgPIE formats of the transitions route need it", {
  expect_error(calcLandReport("downscaledmagpieSplit", "magpie", c(2015, 2050), NULL, "fade", "transitions"),
               "needs downscaling")
})
