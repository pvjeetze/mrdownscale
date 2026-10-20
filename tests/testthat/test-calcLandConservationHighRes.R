# Two clusters of two grid cells each, one year; land conservation of the model by cluster and WDPA
# data by grid cell (base protection for the policy countries and the others).
landTypes <- c("crop", "past", "forestry", "primforest", "secdforest", "urban", "other")
consvMapping <- data.frame(cell = c("1.1.AAA", "2.2.AAA", "3.3.BBB", "4.4.BBB"),
                           cluster = c("AAA.1", "AAA.1", "BBB.2", "BBB.2"))

consvCase <- function() {
  model <- new.magpie(c("AAA.1", "BBB.2"), "y2000",
                      paste(rep(landTypes, each = 2), c("protect", "restore"), sep = "."), fill = 0)
  data <- new.magpie(consvMapping[["cell"]], "y2000", paste0("WDPA.", landTypes), fill = 0,
                     sets = c("x.y.iso", "year", "base.land"))
  return(list(model = model, data = data))
}

# Mocked sources of a case: WDPA base protection for all countries, no future conservation.
consvSources <- function(x, land = NULL, env = parent.frame()) {
  gms <- list(c22_base_protect = "WDPA", c22_base_protect_noselect = "WDPA", c22_protect_scenario = "none",
              c22_protect_scenario_noselect = "none")
  localMagpieSources(list(config = list(gms = gms), clustermap = consvMapping, wdpa = x[["data"]],
                          landConservation = list(landConservation = x[["model"]],
                                                  policyCountries = c("AAA", "BBB"), fader = NULL)),
                     outputs = list("LandTransitionsHighRes.input=magpie" = list(landHr = land)), env = env)
}

test_that("without land, the data are scaled to the cluster values of the model", {
  x <- consvCase()
  # cluster AAA: the model protects 2 Mha natural vegetation (1 by restoration) and 0.4 Mha pasture
  x[["model"]]["AAA.1", , c("primforest.protect", "other.protect", "other.restore")] <- c(0.5, 0.5, 1)
  x[["model"]]["AAA.1", , "past.protect"] <- 0.4
  x[["data"]]["1.1.AAA", , c("WDPA.primforest", "WDPA.past")] <- c(0.5, 0.2)
  x[["data"]]["2.2.AAA", , "WDPA.other"] <- 0.5
  # cluster BBB: data without protection in the model
  x[["data"]]["3.3.BBB", , "WDPA.other"] <- 1
  consvSources(x)
  out <- calcLandConservationHighRes("magpie", onLand = FALSE)[["x"]]
  expect_equal(as.vector(out["1.1.AAA", , c("primforest", "past")]), c(1, 0.4))
  expect_equal(as.vector(out["2.2.AAA", , "other"]), 1)
  expect_equal(sum(out[c("3.3.BBB", "4.4.BBB"), , ]), 0)
})

test_that("on land, protection stays within the land of its type and meets the model where possible", {
  x <- consvCase()
  x[["model"]]["AAA.1", , "primforest.protect"] <- 1.5
  x[["model"]]["AAA.1", , "past.protect"] <- 0.05
  # grid cell 1 has less natural vegetation (1 Mha) than the data protect (1.5 Mha)
  x[["data"]]["1.1.AAA", , c("WDPA.primforest", "WDPA.past")] <- c(1.5, 0.1)
  x[["data"]]["2.2.AAA", , "WDPA.other"] <- 0.2
  land <- new.magpie(consvMapping[["cell"]], "y2000", paste0(landTypes, ".total"), fill = 0)
  land["1.1.AAA", , c("primforest", "secdforest", "past")] <- c(0.6, 0.4, 0.1)
  land["2.2.AAA", , c("other", "past")] <- c(2, 1)
  consvSources(x, land)
  out <- calcLandConservationHighRes("magpie", onLand = TRUE)[["x"]]
  # natural vegetation: 1 Mha in grid cell 1 (all of it, split as its land), the rest of 1.5 in
  # grid cell 2 (0.2 + 0.3 on its unprotected land); pasture: the model protects half of the data
  expect_equal(as.vector(out["1.1.AAA", , c("primforest", "secdforest")]), c(0.6, 0.4), tolerance = 1e-5)
  expect_equal(as.vector(out["2.2.AAA", , "other"]), 0.5, tolerance = 1e-5)
  expect_equal(as.vector(out["1.1.AAA", , "past"]), 0.05, tolerance = 1e-5)
  expect_equal(sum(out[c("3.3.BBB", "4.4.BBB"), , ]), 0)
})

test_that("on land, protection beyond the data fills grid cells with protection in the data first", {
  x <- consvCase()
  land <- new.magpie(consvMapping[["cell"]], "y2000", paste0(landTypes, ".total"), fill = 0)
  # cluster BBB: other land of 1 Mha in both grid cells, data protection of 0.2 Mha in grid cell 3 only
  land[c("3.3.BBB", "4.4.BBB"), , "other"] <- 1
  x[["data"]]["3.3.BBB", , "WDPA.other"] <- 0.2
  # the model protects 1.5 Mha: grid cell 3 is full, the rest goes to grid cell 4
  x[["model"]]["BBB.2", , "other.protect"] <- 1.5
  consvSources(x, land)
  out <- calcLandConservationHighRes("magpie", onLand = TRUE)[["x"]]
  expect_equal(as.vector(out[c("3.3.BBB", "4.4.BBB"), , "other"]), c(1, 0.5), tolerance = 1e-4)
})

test_that("base protection follows the WDPA data by year, the last year after it", {
  x <- consvCase()
  x[["model"]] <- new.magpie(c("AAA.1", "BBB.2"), c("y1995", "y2000", "y2005"),
                             paste(rep(landTypes, each = 2), c("protect", "restore"), sep = "."), fill = 0)
  x[["model"]]["AAA.1", , "primforest.protect"] <- 1
  x[["data"]] <- new.magpie(consvMapping[["cell"]], c("y1995", "y2000"), paste0("WDPA.", landTypes), fill = 0,
                            sets = c("x.y.iso", "year", "base.land"))
  x[["data"]][c("1.1.AAA", "2.2.AAA"), , "WDPA.primforest"] <- c(0.2, 0.8, 0.5, 0.5)
  consvSources(x)
  out <- calcLandConservationHighRes("magpie", onLand = FALSE)[["x"]]
  # grid cell 1 holds 20 % of the cluster data in 1995 and 50 % from 2000 on (taking 2000 for all
  # years would give 50 % in 1995)
  expect_equal(as.vector(out["1.1.AAA", , "primforest"]), c(0.2, 0.5, 0.5))
})

test_that("the data are matched by grid cell", {
  x <- consvCase()
  x[["model"]]["AAA.1", , "primforest.protect"] <- 1
  x[["data"]]["1.1.AAA", , "WDPA.primforest"] <- 0.5
  x[["data"]]["2.2.AAA", , "WDPA.primforest"] <- 1.5
  consvSources(x)
  expected <- calcLandConservationHighRes("magpie", onLand = FALSE)[["x"]]
  # the data in another grid cell order
  x[["data"]] <- x[["data"]][rev(consvMapping[["cell"]]), , ]
  consvSources(x)
  out <- calcLandConservationHighRes("magpie", onLand = FALSE)[["x"]]
  expect_equal(getCells(out), consvMapping[["cell"]])
  expect_equal(as.vector(out), as.vector(expected))
})
