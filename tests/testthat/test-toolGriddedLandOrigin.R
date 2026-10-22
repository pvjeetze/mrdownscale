# transitions (grid cells x years x "stage.from.to"), values grid cell by grid cell, year by year, item by item
stageTransitions <- function(cells, years, names, values) {
  x <- new.magpie(cells, years, names, fill = values)
  getSets(x, fulldim = FALSE)[3] <- "stage.from.to"
  return(x)
}

test_that("outflow takes the newest area first or all area in proportion", {
  initial <- new.magpie("A.1", "y2000", c("crop", "other"), c(2, 1))
  transitions <- stageTransitions("A.1", c("y2010", "y2020"), c("change.crop.other", "change.other.crop"),
                                  c(0.5, 0, 0, 0.7))

  newest <- toolGriddedLandOrigin(transitions, initial, "other")
  expect_equal(getNames(newest), c("other.other", "other.crop"))
  expect_equal(as.vector(newest[, , "other.crop"]), c(0, 0.5, 0))
  expect_equal(as.vector(newest[, , "other.other"]), c(1, 1, 0.8))

  proportional <- toolGriddedLandOrigin(transitions, initial, "other", order = "proportional")
  expect_equal(as.vector(proportional[, "y2020", "other.crop"]), 0.5 * (1 - 0.7 / 1.5))
  expect_equal(as.vector(proportional[, "y2020", "other.other"]), 1 - 0.7 / 1.5)
})

test_that("area keeps its origin between tracked types, with the stages in order", {
  initial <- new.magpie("A.1", "y2000", c("crop", "othernat", "youngsecdf"), c(1, 1, 0))
  transitions <- stageTransitions("A.1", c("y2010", "y2020"),
                                  c("recovery.othernat.youngsecdf", "maturation.youngsecdf.secdforest",
                                    "optimisation.crop.othernat"),
                                  c(0, 0.3, 0, 0.1, 0.4, 0))
  out <- toolGriddedLandOrigin(transitions, initial, c("othernat", "youngsecdf"))
  expect_equal(as.vector(out[, "y2020", c("othernat.othernat", "othernat.crop")]), c(1, 0.1))
  expect_equal(as.vector(out[, "y2020", c("youngsecdf.othernat", "youngsecdf.crop")]), c(0, 0.2))
})

test_that("outflow beyond the start area passes through from the inflow of the stage", {
  initial <- new.magpie("A.1", "y2000", c("urban", "other", "past"), c(1, 0.1, 1))
  transitions <- stageTransitions("A.1", "y2010", c("change.urban.other", "change.other.past"), c(0.3, 0.3))
  out <- expect_silent(toolGriddedLandOrigin(transitions, initial, "other"))
  expect_equal(as.vector(out[, "y2010", c("other.other", "other.urban")]), c(0, 0.1))
})

test_that("the origins sum to the states of the tracked types", {
  set.seed(4)
  types <- c("crop", "past", "other")
  years <- c("y2000", "y2010", "y2020", "y2030")
  states <- new.magpie(c("A.1", "A.2", "B.3"), years, types, fill = runif(3 * 4 * 3))
  # the same area of each grid cell in all years
  area <- setYears(dimSums(states[, 1, ], dim = 3), NULL)
  for (year in years[-1]) {
    states[, year, ] <- states[, year, ] / dimSums(states[, year, ], dim = 3) * area
  }

  # annual transitions of the 10-year steps, as transitions of the step labelled with its end year
  transitions <- setYears(toolMinimalTransitions(states) * 10, years[-1])
  transitions <- add_dimension(transitions, dim = 3.1, add = "stage", nm = "change")
  getSets(transitions, fulldim = FALSE)[3] <- "stage.from.to"
  for (order in c("newest", "proportional")) {
    out <- toolGriddedLandOrigin(transitions, states[, 1, ], c("past", "other"), order = order)
    total <- dimSums(out, dim = "origin")
    for (k in c("past", "other")) {
      expect_equal(as.vector(total[, , k]), as.vector(states[, , k]), tolerance = 1e-12)
    }
  }
})
