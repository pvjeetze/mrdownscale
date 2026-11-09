# A run with one cluster A.1 of two grid cells, three years (y1995 the start) and secondary forest that
# matures at age 5:
# - grid cell 1: y2000 cropland 1 to othernat in the optimisation, half of it restoration (bound 0.5 of the
#   cluster), and 0.2 to secdforest, all of it restoration (bound 0.2; in MAgPIE cropland becomes secondary
#   forest only by restoration); y2005 othernat 0.6 to cropland;
# - grid cell 2 without transitions: cropland 1 (0.5 food crops, 0.3 bioenergy crops, 0.2 fallow), pasture 1,
#   forestry aff 0.3, ndc 0.2, plant 0.5.
btc2Case <- function() {
  cells <- c("1.1.A", "2.2.A")
  years <- c("y1985", "y1995", "y2000", "y2005")
  main <- c("crop", "past", "forestry", "primforest", "secdforest", "urban", "other")
  landHr <- new.magpie(cells, years, main, fill = 0)
  landHr["1.1.A", , "crop"] <- c(2, 2, 0.8, 1.4)
  landHr["1.1.A", , "secdforest"] <- c(1, 1, 1.2, 1.2)
  landHr["1.1.A", , "other"] <- c(1, 1, 2, 1.4)
  landHr["2.2.A", , c("crop", "past", "forestry")] <- 1

  split <- c("crop.kfo.rainfed", "crop.kbe.rainfed", "crop.fallow.total", "past.total.total", "forestry.aff.total",
             "forestry.ndc.total", "forestry.plant.total", "primforest.total.total", "secdforest.total.total",
             "urban.total.total", "other.othernat.total", "other.youngsecdf.total")
  splitHr <- new.magpie(cells, years, split, fill = 0, sets = c("x", "y", "iso", "year", "land", "type", "water"))
  splitHr["1.1.A", , "crop.kfo.rainfed"] <- landHr["1.1.A", , "crop"]
  splitHr["1.1.A", , "secdforest.total.total"] <- landHr["1.1.A", , "secdforest"]
  splitHr["1.1.A", , "other.othernat.total"] <- landHr["1.1.A", , "other"]
  splitHr["2.2.A", , ] <- c(0.5, 0.3, 0.2, 1, 0.3, 0.2, 0.5, 0, 0, 0, 0, 0)[rep(1:12, each = 4)]

  transitions <- c("optimisation.crop.total.other.othernat", "optimisation.crop.total.secdforest.total",
                   "optimisation.other.othernat.crop.total")
  transitionsHr <- new.magpie(cells, years[-1], transitions, fill = 0,
                              sets = c("x", "y", "iso", "year", "stage", "from", "fromType", "to", "toType"))
  transitionsHr["1.1.A", "y2000", transitions[1:2]] <- c(1, 0.2)
  transitionsHr["1.1.A", "y2005", transitions[3]] <- 0.6

  restorationLr <- new.magpie("A.1", years[-1], c("secdforest", "other"), fill = 0)
  restorationLr[, "y2000", ] <- c(0.2, 0.5)
  densityLr <- new.magpie("A.1", years[-1], c("ac0", "ac5", "ac10"), fill = c(0, 0, 0, 25, 25, 25, 40, 40, 40))
  return(list(landHr = landHr, splitHr = splitHr, transitionsHr = transitionsHr, restorationLr = restorationLr,
              densityLr = densityLr, mapping = data.frame(cell = cells, cluster = "A.1")))
}

# A variant of the case without restoration (all transitions are unrestored), years y1985, y1995, y2000, y2005:
# - grid cell 1: y2000 cropland 0.4 to secondary forest (secdforest 1 to 1.4); y2005 secondary forest 0.7 to
#   cropland;
# - grid cell 2: y2000 pasture 0.3, plant 0.2 and other_planted 0.1 to othernat; y2005 pasture 0.2 to othernat
#   and othernat 0.4 to cropland in the same stage.
btc2OriginCase <- function() {
  case <- btc2Case()
  case[["restorationLr"]][, , ] <- 0
  landHr <- case[["landHr"]]
  landHr[, , ] <- 0
  landHr["1.1.A", , c("crop", "secdforest", "other")] <- c(2, 2, 1.6, 2.3, 1, 1, 1.4, 0.7, 1, 1, 1, 1)
  landHr["2.2.A", , c("crop", "past", "forestry", "other")] <-
    c(1, 1, 1, 1.4, 1, 1, 0.7, 0.5, 1, 1, 0.7, 0.7, 0, 0, 0.6, 0.4)
  splitHr <- mbind(case[["splitHr"]], setNames(case[["splitHr"]][, , "forestry.plant.total"] * 0,
                                               "forestry.other_planted.total"))
  splitHr[, , ] <- 0
  splitHr["1.1.A", , "crop.kfo.rainfed"] <- landHr["1.1.A", , "crop"]
  splitHr["1.1.A", , "secdforest.total.total"] <- landHr["1.1.A", , "secdforest"]
  splitHr["1.1.A", , "other.othernat.total"] <- landHr["1.1.A", , "other"]
  splitHr["2.2.A", , c("crop.kfo.rainfed", "crop.kbe.rainfed", "crop.fallow.total")] <-
    c(0.5, 0.5, 0.5, 0.9, rep(0.3, 4), rep(0.2, 4))
  splitHr["2.2.A", , c("past.total.total", "forestry.aff.total", "forestry.ndc.total")] <-
    c(1, 1, 0.7, 0.5, rep(0.3, 4), rep(0.2, 4))
  splitHr["2.2.A", , c("forestry.plant.total", "forestry.other_planted.total", "other.othernat.total")] <-
    c(0.4, 0.4, 0.2, 0.2, 0.1, 0.1, 0, 0, 0, 0, 0.6, 0.4)

  transitions <- c("optimisation.crop.total.secdforest.total", "optimisation.secdforest.total.crop.total",
                   "optimisation.past.total.other.othernat", "optimisation.forestry.plant.other.othernat",
                   "optimisation.forestry.other_planted.other.othernat", "optimisation.other.othernat.crop.total")
  transitionsHr <- new.magpie(getCells(landHr), getYears(landHr)[-1], transitions, fill = 0,
                              sets = c("x", "y", "iso", "year", "stage", "from", "fromType", "to", "toType"))
  transitionsHr["1.1.A", "y2000", transitions[1]] <- 0.4
  transitionsHr["1.1.A", "y2005", transitions[2]] <- 0.7
  transitionsHr["2.2.A", , transitions[3]] <- c(0, 0.3, 0.2)
  transitionsHr["2.2.A", "y2000", transitions[4:5]] <- c(0.2, 0.1)
  transitionsHr["2.2.A", "y2005", transitions[6]] <- 0.4
  case[["landHr"]] <- landHr
  case[["splitHr"]] <- splitHr
  case[["transitionsHr"]] <- transitionsHr
  return(case)
}

# The case as synthetic sources and outputs of calcLandReportBtC2.
localBtC2Sources <- function(case, affPlantation = 0, env = parent.frame()) {
  .reportKey <- function(outputFormat) {
    return(paste(c("LandReport", paste0("outputFormat=", outputFormat), "input=magpie",
                   "harmonizationPeriod=2015,2050", paste0("yearsSubset=", paste(seq(1995, 2100, 5), collapse = ",")),
                   "harmonization=none", "downscaling=transitions"), collapse = "."))
  }
  outputs <- list(list(transitionsHr = case[["transitionsHr"]]), case[["landHr"]], case[["splitHr"]])
  names(outputs) <- c("LandTransitionsHighRes.input=magpie", .reportKey("downscaledmagpie"),
                      .reportKey("downscaledmagpieSplit"))
  localMagpieSources(list(clustermap = case[["mapping"]], config = list(gms = list(s32_aff_plantation = affPlantation)),
                          restoration = case[["restorationLr"]], secdforestCarbonDensity = case[["densityLr"]]),
                     outputs = outputs, env = env)
}

test_that("the state classes map the land split, afforestation by s32_aff_plantation", {
  localBtC2Sources(btc2Case())
  out <- calcLandReportBtC2()[["x"]]
  expect_equal(getNames(out), c("cropland_other", "cropland_2Gbioen", "grassland", "nat_regen_forest",
                                "forest_planted", "rest_abn_crop", "rest_abn_grass", "rest_forest", "other",
                                "built_up_areas"))
  expect_equal(getYears(out), c("y1995", "y2000", "y2005"))
  expect_equal(as.vector(out["2.2.A", "y2005", c("cropland_other", "cropland_2Gbioen", "grassland",
                                                 "nat_regen_forest", "forest_planted")]), c(0.7, 0.3, 1, 0.5, 0.5))
  localBtC2Sources(btc2Case(), affPlantation = 1)
  out <- calcLandReportBtC2()[["x"]]
  expect_equal(as.vector(out["2.2.A", "y2005", c("nat_regen_forest", "forest_planted")]), c(0.2, 0.8))
})

test_that("abandoned and restored cropland are rest_abn_crop until conversion or maturation", {
  localBtC2Sources(btc2Case())
  out <- calcLandReportBtC2()[["x"]]
  # y2000: in other land 0.5 abandoned and 0.5 restored cropland; in secondary forest 0.2 entered directly
  expect_equal(as.vector(out["1.1.A", "y2000", c("rest_abn_crop", "other", "nat_regen_forest")]), c(1.2, 1, 1))
  # y2005: the conversion takes the abandoned cropland and 0.1 of the initial other land, not the restored
  # cropland; the direct secondary forest has matured (age 5)
  expect_equal(as.vector(out["1.1.A", "y2005", c("rest_abn_crop", "other", "nat_regen_forest")]), c(0.5, 0.9, 1.2))
  expect_equal(as.vector(dimSums(out, dim = 3)), as.vector(dimSums(btc2Case()[["landHr"]][, -1, ], dim = 3)))
})

test_that("restoration is marked in all stages of the model's transitions", {
  # the transition of cropland to secondary forest in the SNV relocation stage (magpie4::landTransitions)
  case <- btc2Case()
  items <- getNames(case[["transitionsHr"]])
  items[items == "optimisation.crop.total.secdforest.total"] <- "snvRelocation.crop.total.secdforest.total"
  getNames(case[["transitionsHr"]]) <- items
  localBtC2Sources(case)
  expect_no_message(out <- calcLandReportBtC2()[["x"]], message = "Check failed")
  expect_equal(as.vector(out["1.1.A", "y2000", c("rest_abn_crop", "nat_regen_forest")]), c(1.2, 1))
})

test_that("conversions take restored land only beyond all other area, by the model as a note", {
  # the cluster has no other unrestored land: the cluster sums take the restored land as well
  case <- btc2Case()
  case[["transitionsHr"]]["1.1.A", "y2005", "optimisation.other.othernat.crop.total"] <- 1.8
  case[["landHr"]]["1.1.A", "y2005", c("crop", "other")] <- c(2.6, 0.2)
  case[["splitHr"]]["1.1.A", "y2005", c("crop.kfo.rainfed", "other.othernat.total")] <- c(2.6, 0.2)
  localBtC2Sources(case)
  expect_no_warning(expect_message(out <- calcLandReportBtC2()[["x"]],
                                   "MAgPIE converted restored land .* up to 0.3 Mha in a time step \\(y2005\\)"))
  expect_equal(as.vector(out["1.1.A", "y2005", c("rest_abn_crop", "other")]), c(0.2, 0))
})

test_that("restored land taken only by the grid cell placement gives a warning", {
  # grid cell 2 holds other land 1 throughout: the cluster sums take no restored land
  case <- btc2Case()
  case[["transitionsHr"]]["1.1.A", "y2005", "optimisation.other.othernat.crop.total"] <- 1.8
  case[["landHr"]]["1.1.A", "y2005", c("crop", "other")] <- c(2.6, 0.2)
  case[["splitHr"]]["1.1.A", "y2005", c("crop.kfo.rainfed", "other.othernat.total")] <- c(2.6, 0.2)
  case[["landHr"]]["2.2.A", , "other"] <- 1
  case[["splitHr"]]["2.2.A", , "other.othernat.total"] <- 1
  localBtC2Sources(case)
  expect_warning(out <- calcLandReportBtC2()[["x"]],
                 "placement converted restored land beyond the model, up to 0.3 Mha in a time step \\(y2005\\)")
  expect_equal(as.vector(out["1.1.A", "y2005", c("rest_abn_crop", "other")]), c(0.2, 0))
  expect_equal(as.vector(out["2.2.A", "y2005", "other"]), 1)
})

test_that("direct secondary forest, origins from pasture and timber plantations, and no maturation below 20 tC/ha", {
  case <- btc2OriginCase()
  classes <- c("rest_abn_crop", "rest_abn_grass", "rest_forest", "other", "nat_regen_forest")
  # maturation at age 5 (carbon density 25 tC/ha in ac5): the direct secondary forest of y2000 has matured in y2005
  localBtC2Sources(case)
  out <- calcLandReportBtC2()[["x"]]
  # y2005, grid cell 1: nothing recovering, other 1, secondary forest 0.7 (the loss of 0.7 from 1.4 takes nothing
  # from the matured land)
  expect_equal(as.vector(out["1.1.A", "y2005", classes[c(1, 4, 5)]]), c(0, 1, 0.7))

  # carbon density never above 20 tC/ha: the direct secondary forest never matures
  case[["densityLr"]][, , ] <- c(0, 0, 0, 15, 15, 15, 20, 20, 20)
  localBtC2Sources(case)
  out <- calcLandReportBtC2()[["x"]]
  # (a) y2000 grid cell 1: 0.4 entered secondary forest from cropland, secondary forest 1.4, so nat_regen_forest
  # 1.4 - 0.4 = 1 and other 1 + 0.4 - 0.4 = 1
  expect_equal(as.vector(out["1.1.A", "y2000", classes[c(1, 4, 5)]]), c(0.4, 1, 1))
  # y2005: the loss of 0.7 of 1.4 takes half of the direct area, 0.4 * (1 - 0.7 / 1.4) = 0.2, still recovering
  # land (c); nat_regen_forest 0.7 - 0.2 = 0.5, other 1 + 0.2 - 0.2 = 1
  expect_equal(as.vector(out["1.1.A", "y2005", classes[c(1, 4, 5)]]), c(0.2, 1, 0.5))
  expect_equal(as.vector(out["1.1.A", "y2005", "cropland_other"]), 2.3)

  # (b) y2000 grid cell 2: pasture 0.3 and plant 0.2 + other_planted 0.1 = 0.3 timber entered othernat
  expect_equal(as.vector(out["2.2.A", "y2000", classes[1:4]]), c(0, 0.3, 0.3, 0))
  expect_equal(as.vector(out["2.2.A", "y2000", c("grassland", "forest_planted", "nat_regen_forest")]),
               c(0.7, 0.2, 0.5))
  # y2005: the outflow of 0.4 takes the area at the start of the stage (0.6, the entries of y2000, 2/3 of each
  # origin), so pasture 0.3 * (1 / 3) = 0.1 and timber 0.3 * (1 / 3) = 0.1 remain; the inflow of 0.2 pasture
  # is added: pasture 0.1 + 0.2 = 0.3, timber 0.1; other land 0.4 is all recovering
  expect_equal(as.vector(out["2.2.A", "y2005", classes[1:4]]), c(0, 0.3, 0.1, 0))
  expect_equal(as.vector(out["2.2.A", "y2005", c("grassland", "cropland_other")]), c(0.5, 1.1))
})
