test_that("cluster states, transitions by stage and sub-type, and output names follow the run", {
  case <- magpieCase()
  localMagpieSources(magpieSources(case))
  out <- calcLandTransitionsHighRes("magpie")[["x"]]
  expect_lt(out[["deviationLr"]], 1e-6)
  expect_equal(getYears(out[["landHr"]]), c("y1985", "y2000", "y2005"))
  expect_equal(getNames(out[["landHr"]]),
               c("crop.total", "past.total", "forestry.aff", "forestry.plant", "primforest.total",
                 "secdforest.total", "urban.total", "other.othernat", "other.youngsecdf"))
  # cell area stays the same
  area <- dimSums(out[["landHr"]], dim = 3)
  expect_equal(as.vector(area[, "y2005", ]), as.vector(area[, "y1985", ]))

  tr <- out[["transitionsHr"]]
  expect_equal(getSets(tr, fulldim = FALSE)[3], "stage.from.fromType.to.toType")
  expect_equal(clusterFlow(tr, "y2000", "disturbance.primforest.total.secdforest.total"), 0.1)
  expect_equal(clusterFlow(tr, "y2000", "maturation.other.youngsecdf.secdforest.total"), 0.05)
  expect_equal(clusterFlow(tr, "y2000", "optimisation.past.total.crop.total"), 0.3)
  # sub-types from the module flows: othernat to crop, crop to plant, aff to othernat, harvest of youngsecdf
  expect_equal(clusterFlow(tr, "y2000", "optimisation.other.othernat.crop.total"), 0.2, tolerance = 1e-6)
  expect_equal(clusterFlow(tr, "y2000", "optimisation.crop.total.forestry.plant"), 0.1, tolerance = 1e-6)
  expect_equal(clusterFlow(tr, "y2000", "optimisation.forestry.aff.other.othernat"), 0.05, tolerance = 1e-6)
  expect_equal(clusterFlow(tr, "y2000", "optimisation.other.youngsecdf.other.othernat"), 0.02, tolerance = 1e-6)
  # the main transitions are the sum over stages and sub-types
  mainTr <- dimSums(tr, dim = c("stage", "fromType", "toType"))
  expect_equal(sum(mainTr[, "y2005", "crop.past"]), 0.1)
})

test_that("the initial state is split into sub-types by the cluster shares", {
  case <- magpieCase()
  localMagpieSources(magpieSources(case))
  land <- calcLandTransitionsHighRes("magpie")[["x"]][["landHr"]]
  # cluster start: aff 0.2 and plant 0.3 of forestry 0.5; youngsecdf 0.2 of other 3
  expect_equal(as.vector(land[, "y1985", "forestry.aff"]), c(0.2, 0, 0.3) * 0.4, tolerance = 1e-6)
  expect_equal(as.vector(land[, "y1985", "other.youngsecdf"]), c(0.7, 0.8, 1.5) * 0.2 / 3, tolerance = 1e-6)
})

test_that("recovery is placed on the othernat inflow from forestry of the previous year", {
  case <- magpieCase()
  localMagpieSources(magpieSources(case))
  tr <- calcLandTransitionsHighRes("magpie")[["x"]][["transitionsHr"]]
  recovery <- as.vector(tr[, "y2005", "recoveryForestry.other.othernat.other.youngsecdf"])
  expect_equal(sum(recovery), 0.04)
  # cell 2 has no forestry, so it had no abandonment
  expect_lt(recovery[2], 1e-4)
  inflow <- as.vector(tr[, "y2000", "optimisation.forestry.aff.other.othernat"])
  expect_true(all(inflow[c(1, 3)] > 0))
})

test_that("cropland stays within the available cropland of the grid cells", {
  case <- magpieCase()
  # cluster cropland 2.1 and 2.0 Mha; grid cell 2 can hold at most 0.2 Mha
  localMagpieSources(magpieSources(case, availableCropland = c(1.5, 0.2, 1)))
  out <- calcLandTransitionsHighRes("magpie")[["x"]]
  crop <- out[["landHr"]][, c("y2000", "y2005"), "crop"]
  # the caps hold to the convergence tolerance of 1e-6 Mha per cluster
  expect_true(all(as.vector(crop) <= rep(c(1.5, 0.2, 1), 2) + 1e-5))
  expect_lt(out[["deviationLr"]], 1e-6)
})

test_that("land reduction on drained peat follows the drained peat, and peat equals toolDownscalePeatland", {
  case <- magpieCase()
  cells <- case[["mapping"]][["cell"]]
  # the drained peat under cropland of the cluster declines by 0.1 in y2005
  peatlandLr <- new.magpie("A.1", c("y2000", "y2005"), c("crop", "intact", "rewetted"),
                           fill = c(0.4, 0.3, 0.5, 0.5, 0, 0.1))
  peatlandStartHr <- new.magpie(cells, "y1995", c("crop", "intact"), fill = c(0.4, 0, 0, 0.5, 0, 0))
  case <- addStage(case, "peatRewetting.crop.past", "y2005", 0.1)
  sources <- magpieSources(case, config = list(peatland = "v2"))
  sources[["MagpieFulldataGdx.peatland"]] <- peatlandLr
  sources[["MagpieRunFiles.peatland"]] <- peatlandStartHr
  localMagpieSources(sources)
  out <- calcLandTransitionsHighRes("magpie")[["x"]]
  rewetting <- as.vector(out[["transitionsHr"]][, "y2005", "peatRewetting.crop.total.past.total"])
  expect_equal(rewetting, c(0.1, 0, 0), tolerance = 1e-5)
  peatHr <- out[["peatlandHr"]]
  expect_equal(as.vector(peatHr[1, "y2005", ]), c(0.3, 0.5, 0.1), tolerance = 1e-5)
  expect_lt(out[["deviationLr"]], 1e-6)

  rewettingHr <- dimSums(out[["transitionsHr"]][, , "peatRewetting"], dim = c("stage", "fromType", "to", "toType"))
  standalone <- toolDownscalePeatland(out[["landHr"]][, c("y2000", "y2005"), ], peatlandLr, peatlandStartHr,
                                      case[["mapping"]], rewettingHr)
  expect_equal(as.vector(standalone[["peatlandHr"]]), as.vector(peatHr))
})

test_that("the SNV relocation follows the relocation target of the grid cells, matched by grid cell", {
  case <- magpieCase()
  # y2005: 0.1 of cropland relocated to othernat (q29_land_snv_trans)
  case[["dataLr"]][["land"]][, "y2005", c("crop", "other")] <- c(1.9, 2.9)
  case[["dataLr"]][["other"]][, "y2005", "othernat"] <- 2.63
  case <- addStage(case, "snvRelocation.crop.othernat", "y2005", 0.1)
  sources <- magpieSources(case)
  # SNV share 0.2 (the first target): the target lies in grid cell 3 only, given in reverse grid cell order
  sources[["snv"]][["shareIso"]][["A"]] <- 0.2
  sources[["snvTarget"]] <- new.magpie(rev(case[["mapping"]][["cell"]]), NULL,
                                       c("SNV20TargetCropland", "SNV50TargetCropland"), fill = c(0.5, 0, 0, 0.5, 0, 0))
  localMagpieSources(sources)
  out <- calcLandTransitionsHighRes("magpie")[["x"]]
  relocation <- as.vector(out[["transitionsHr"]][, "y2005", "snvRelocation.crop.total.other.othernat"])
  expect_equal(sum(relocation), 0.1, tolerance = 1e-6)
  expect_gt(relocation[3], 0.1 - 1e-4)
  expect_lt(out[["deviationLr"]], 1e-6)
})

test_that("cropland is split by the cluster shares; a cluster without cropland takes the shares of all", {
  case <- twoClusters(magpieCase())
  # cluster B has no crop composition in y2005
  case[["cropSplitLr"]]["B.2", "y2005", ] <- 0
  case[["cropAreaLr"]]["B.2", "y2005", ] <- 0
  localMagpieSources(magpieSources(case))
  out <- calcLandTransitionsHighRes("magpie")[["x"]]
  cropHr <- dimSums(out[["landHr"]][, , "crop"], dim = 3)
  split <- out[["splitHr"]][["crop"]]
  area <- out[["splitHr"]][["crop_area"]]
  expect_equal(getYears(split), c("y1985", "y2000", "y2005"))
  # every grid cell gets the composition of its cluster; the start that of the first year
  expect_equal(as.vector(split["1.1.A", "y1985", ]), as.vector(cropHr["1.1.A", "y1985", ]) * c(6, 2, 1) / 9)
  expect_equal(as.vector(area["2.2.A", "y2000", ]),
               as.vector(cropHr["2.2.A", "y2000", ]) * 6 / 9 * c(0.75, 0.25))
  # cluster B in y2005: the composition of all clusters, that is of cluster A
  expect_equal(as.vector(split["4.4.B", "y2005", ]), as.vector(cropHr["4.4.B", "y2005", ]) * c(6, 2, 1) / 9)
  expect_equal(as.vector(dimSums(area, dim = 3)), as.vector(split[, , "crop_area"]))
})

test_that("other harmonisations stop with a message", {
  case <- magpieCase()
  localMagpieSources(magpieSources(case))
  expect_error(calcLandTransitionsHighRes("magpie", harmonizationPeriod = c(2015, 2050), harmonization = "fade"),
               "only \"none\"")
})
