test_that("the operations follow the model order, with disturbance only where the model has it", {
  case <- ageCase()
  localMagpieSources(case[["sources"]])
  x <- calcAgeClassesInput("magpie")[["x"]]
  expect_equal(x[["pools"]], c("secdforest.total", "forestry.plant", "other.othernat", "other.youngsecdf"))
  expect_equal(x[["groups"]][["forestry"]], "forestry.plant")
  kinds <- vapply(x[["operations"]], function(o) o[["kind"]], character(1))
  # no disturbance of forestry in the case: within (secondary forest), transfer (recovery), shift, transfer
  # (maturation), one exchange per land type
  expect_equal(kinds, c("within", "transfer", "shift", "transfer", rep("exchange", 4)))
  expect_equal(x[["operations"]][[2]][["stages"]], c("recoveryForestry", "recovery"))
  expect_equal(getItems(x[["operations"]][[5]][["outflowLr"]], dim = 3), c("ac0", "ac5", "ac10", "acx"))
  expect_equal(as.vector(x[["operations"]][[5]][["outflowLr"]][, "y2005", "acx"]), 0.5)
})

test_that("forestry types with disturbance get a within operation", {
  case <- ageCase()
  flows <- case[["sources"]][["ageClassTransitions"]]
  disturbance <- flows[, , "disturbance.secdforest.total"]
  getItems(disturbance, dim = 3.3, raw = TRUE) <- "plant"
  getItems(disturbance, dim = 3.2, raw = TRUE) <- "forestry"
  case[["sources"]][["ageClassTransitions"]] <- mbind(flows, disturbance)
  localMagpieSources(case[["sources"]])
  operations <- calcAgeClassesInput("magpie")[["x"]][["operations"]]
  expect_equal(operations[[2]][["kind"]], "within")
  expect_equal(operations[[2]][["pool"]], "forestry.plant")
})
