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
# The switches of module 22 can be changed by gms, the clustermap by mapping.
consvSources <- function(x, land = NULL, env = parent.frame(), gms = list(), mapping = consvMapping) {
  gmsAll <- list(c22_base_protect = "WDPA", c22_base_protect_noselect = "WDPA", c22_protect_scenario = "none",
                 c22_protect_scenario_noselect = "none")
  gmsAll[names(gms)] <- gms
  localMagpieSources(list(config = list(gms = gmsAll), clustermap = mapping, wdpa = x[["data"]],
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

test_that("grid cells of countries outside the policy follow the noselect base protection", {
  # no base protection in the policy countries (AAA, BBB), WDPA data in all others
  x <- consvCase()
  mapping <- rbind(consvMapping, data.frame(cell = "5.5.CCC", cluster = "CCC.3"))
  x[["model"]] <- mbind(x[["model"]], new.magpie("CCC.3", "y2000", getNames(x[["model"]]), fill = 0))
  x[["data"]] <- mbind(x[["data"]], new.magpie("5.5.CCC", "y2000", getNames(x[["data"]]), fill = 0,
                                               sets = c("x.y.iso", "year", "base.land")))
  # the model protects the 0.5 Mha primary forest of the data in grid cell 5 (country CCC) and in grid cell 1 (AAA)
  x[["model"]][c("AAA.1", "CCC.3"), , "primforest.protect"] <- 0.5
  x[["data"]][c("1.1.AAA", "5.5.CCC"), , "WDPA.primforest"] <- 0.5
  consvSources(x, gms = list(c22_base_protect = "none"), mapping = mapping)
  out <- calcLandConservationHighRes("magpie", onLand = FALSE)[["x"]]
  # grid cell 5: base protection of the data (noselect), scaling factor 0.5 / 0.5 = 1, so 0.5 Mha
  expect_equal(as.vector(out["5.5.CCC", , "primforest"]), 0.5)
  # grid cell 1 (policy country): base protection none, so no data and 0 Mha
  expect_equal(as.vector(out["1.1.AAA", , "primforest"]), 0)
})
