# Two grid cells on the 0.5 degree grid and two years, land in Mha by class. Grid cell A at lon -179.75, lat 89.75
# (first index of both axes), grid cell B at lon 10.25, lat -20.75.
# - A: area 2.5 Mha; y2000 0.5, 0.75, 1.25; y2005 1, 0.5, 1
# - B: area 0.8 Mha; y2000 0.2, 0.2, 0.4; y2005 0.4, 0.2, 0.2
btc2Classes <- c("cropland_other", "grassland", "other")
btc2Land <- function() {
  landHr <- new.magpie(c("-179p75.89p75.AAA", "10p25.-20p75.BBB"), c("y2000", "y2005"), btc2Classes,
                       fill = c(0.5, 0.2, 1, 0.4, 0.75, 0.2, 0.5, 0.2, 1.25, 0.4, 1, 0.2),
                       sets = c("x.y.iso", "year", "class"))
  return(landHr)
}

# Sources of fullBTC2: the land report and the config of a run with the given title.
localFullBtC2Sources <- function(landHr, title, env = parent.frame()) {
  localMagpieSources(list(config = list(title = title)), outputs = list(LandReportBtC2 = landHr), env = env)
  local_mocked_bindings(toolWriteMadratLog = function(...) invisible(NULL), .env = env)
}

test_that("fullBTC2 names the file by the title and the project name", {
  withr::local_dir(withr::local_tempdir())
  landHr <- btc2Land()
  localFullBtC2Sources(landHr, "run1")
  fullBTC2()
  expect_true(file.exists("BtC2v04_run1.nc"))
  # the title contains the project name (in any case)
  localFullBtC2Sources(landHr, "btc2V04_run2")
  fullBTC2()
  expect_true(file.exists("btc2V04_run2.nc"))
  expect_false(file.exists("BtC2v04_btc2V04_run2.nc"))
  localFullBtC2Sources(landHr, "run3")
  fullBTC2(projectName = "other")
  expect_true(file.exists("other_run3.nc"))
})

test_that("fullBTC2 writes the shares of the classes and the grid cell area on the 0.5 degree grid", {
  withr::local_dir(withr::local_tempdir())
  landHr <- btc2Land()
  localFullBtC2Sources(landHr, "run1")
  fullBTC2()
  nc <- ncdf4::nc_open("BtC2v04_run1.nc")
  withr::defer(ncdf4::nc_close(nc))

  # axes: 720 longitudes, 360 latitudes from north to south, the classes, the years
  expect_equal(nc[["dim"]][["lon"]][["len"]], 720)
  expect_equal(nc[["dim"]][["lat"]][["len"]], 360)
  expect_equal(as.vector(ncdf4::ncvar_get(nc, "lon"))[c(1, 720)], c(-179.75, 179.75))
  expect_equal(as.vector(ncdf4::ncvar_get(nc, "lat"))[c(1, 360)], c(89.75, -89.75))
  expect_equal(as.vector(ncdf4::ncvar_get(nc, "lc_class")), 1:3)
  expect_equal(nc[["dim"]][["lc_class"]][["units"]], "1=cropland_other/2=grassland/3=other")
  expect_equal(as.vector(ncdf4::ncvar_get(nc, "time")), c(2000, 2005))

  # grid cell A is at index [1, 1], grid cell B at [(10.25 + 179.75) / 0.5 + 1, (89.75 + 20.75) / 0.5 + 1]
  indexA <- c(1, 1)
  indexB <- c(381, 222)
  share <- ncdf4::ncvar_get(nc, "LC_area_share")
  expect_equal(dim(share), c(720, 360, 3, 2))
  # shares: class area / grid cell area, A: 0.5 / 2.5 = 0.2, 0.75 / 2.5 = 0.3, 1.25 / 2.5 = 0.5 (y2000) and
  # 1 / 2.5 = 0.4, 0.5 / 2.5 = 0.2, 1 / 2.5 = 0.4 (y2005); B: 0.2 / 0.8 = 0.25, 0.25, 0.4 / 0.8 = 0.5 (y2000) and
  # 0.4 / 0.8 = 0.5, 0.25, 0.25 (y2005)
  expect_equal(share[indexA[1], indexA[2], , 1], c(0.2, 0.3, 0.5))
  expect_equal(share[indexA[1], indexA[2], , 2], c(0.4, 0.2, 0.4))
  expect_equal(share[indexB[1], indexB[2], , 1], c(0.25, 0.25, 0.5))
  expect_equal(share[indexB[1], indexB[2], , 2], c(0.5, 0.25, 0.25))
  # all other grid cells have no value
  expect_equal(sum(!is.na(share)), 2 * 3 * 2)

  # grid cell area in Mha: 0.5 + 0.75 + 1.25 = 2.5 (A), 0.2 + 0.2 + 0.4 = 0.8 (B)
  area <- ncdf4::ncvar_get(nc, "pixel_area")
  expect_equal(dim(area), c(720, 360))
  expect_equal(c(area[indexA[1], indexA[2]], area[indexB[1], indexB[2]]), c(2.5, 0.8))
  expect_equal(sum(!is.na(area)), 2)
})

test_that("fullBTC2 writes the land of the grid cells as cell.land_btc2_0.5.mz", {
  withr::local_dir(withr::local_tempdir())
  landHr <- btc2Land()
  localFullBtC2Sources(landHr, "run1")
  fullBTC2()
  expect_true(file.exists("cell.land_btc2_0.5.mz"))
  written <- read.magpie("cell.land_btc2_0.5.mz")
  expect_equal(getYears(written), c("y2000", "y2005"))
  expect_equal(getNames(written), btc2Classes)
  # the values of the input: A y2000 0.5, 0.75, 1.25 and B y2005 0.4, 0.2, 0.2
  expect_equal(as.vector(written[1, "y2000", ]), c(0.5, 0.75, 1.25), tolerance = 1e-5)
  expect_equal(as.vector(written[2, "y2005", ]), c(0.4, 0.2, 0.2), tolerance = 1e-5)
})
