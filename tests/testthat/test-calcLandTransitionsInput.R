test_that("the stages of each year are matrices in the order of the model, the states by sub-type", {
  case <- magpieCase()
  localMagpieSources(magpieSources(case))
  x <- calcLandTransitionsInput("magpie")[["x"]]
  expect_equal(names(x[["transitionsLr"]]), c("y2000", "y2005"))
  expect_equal(names(x[["transitionsLr"]][["y2000"]]),
               c("disturbance", "recoveryForestry", "recovery", "maturation", "optimisation"))
  optimisation <- x[["transitionsLr"]][["y2000"]][["optimisation"]]
  expect_equal(dimnames(optimisation),
               list("A.1", c("past.crop", "othernat.crop", "crop.plant", "aff.othernat", "youngsecdf.othernat",
                             "crop.past")))
  expect_equal(as.vector(optimisation), c(0.3, 0.2, 0.1, 0.05, 0.02, 0))
  expect_equal(as.vector(x[["transitionsLr"]][["y2005"]][["recoveryForestry"]]), 0.04)
  expect_equal(unname(x[["mainOf"]]), c("crop", "past", "forestry", "forestry", "primforest", "secdforest", "urban",
                                        "other", "other"))
  expect_equal(as.vector(x[["statesLr"]][, c("othernat", "youngsecdf"), "y2005"]), c(2.53, 0.27))
  # the start map with forestry and other land in equal parts of their sub-types
  expect_equal(unname(x[["startHr"]][, "aff"]), c(0.2, 0, 0.3) / 2)
  expect_equal(unname(x[["startHr"]][, "past"]), c(0.5, 1, 0.5) + c(0, 0.5, 0))
  expect_equal(names(x[["splitLr"]]), c("crop", "crop_area"))
})

test_that("rules place recovery, and where the run has them the SNV relocation, peat and urban land", {
  case <- magpieCase()
  localMagpieSources(magpieSources(case))
  rules <- calcLandTransitionsInput("magpie")[["x"]][["rules"]]
  expect_equal(vapply(rules, function(r) r[["kind"]], character(1)), c("previousInflow", "previousInflow"))
  expect_equal(rules[[1]][["fromTypes"]], c("aff", "plant"))

  case <- addStage(case, "snvRelocation.crop.othernat", "y2005", 0.1)
  sources <- magpieSources(case, config = list(peatland = "v2", urban = "exo_nov21", c34_urban_scenario = "SSP2"))
  cells <- case[["mapping"]][["cell"]]
  sources[["MagpieFulldataGdx.peatland"]] <- new.magpie("A.1", c("y2000", "y2005"), c("crop", "intact"), fill = 0.1)
  sources[["MagpieRunFiles.peatland"]] <- new.magpie(cells, "y1995", c("crop", "intact"), fill = 0.1)
  sources[["urban"]] <- new.magpie(cells, c("y2000", "y2005"), "SSP2", fill = 0.1)
  # SNV share 0.35, between the 20 % and 50 % targets: the 20 % target plus half of the step to the 50 % target
  sources[["snv"]][["shareIso"]][["A"]] <- 0.35
  sources[["snvTarget"]] <- new.magpie(cells, NULL, c("SNV20TargetCropland", "SNV50TargetCropland"),
                                       fill = c(0.2, 0, 0, 0.4, 0.2, 0))
  localMagpieSources(sources)
  x <- calcLandTransitionsInput("magpie")[["x"]]
  kinds <- stats::setNames(vapply(x[["rules"]], function(r) r[["kind"]], character(1)),
                           vapply(x[["rules"]], function(r) r[["stage"]], character(1)))
  expect_equal(kinds, c(recoveryForestry = "previousInflow", recovery = "previousInflow",
                        snvRelocation = "depletingTarget", peatRewetting = "peatland", optimisation = "targetArea"))
  expect_equal(unname(x[["rules"]][[3]][["targetHr"]]), c(0.3, 0.1, 0))
  expect_equal(x[["rules"]][[4]][["onTypes"]], "crop")
  # no new inflow to secondary forest, youngsecdf and primary forest in the urban compensation
  allowed <- x[["rules"]][[5]][["allowedTransitions"]]
  expect_false(any(allowed[, c("secdforest", "youngsecdf", "primforest")]))
})
