# Synthetic MAgPIE runs for the tests of the downscaling along the transitions.

# A MAgPIE run with one cluster of three cells and two years, consistent with its land accounting, and its
# transitions by stage and sub-type.
# y2000: disturbance 0.1, recovery 0.1, maturation 0.05; past to crop 0.3, othernat to crop 0.2,
# crop to plant 0.1, aff to othernat 0.05, harvest of youngsecdf 0.02 (regrowing as othernat).
# y2005: recovery 0.04 (all from the forestry abandonment of y2000); crop to past 0.1.
# Cropland: crop area 2/3, fallow 2/9, tree cover 1/9; crop area 3/4 rainfed food crops (tece), 1/4
# irrigated bioenergy grass (begr).
magpieCase <- function() {
  cells <- c("1.1.A", "2.2.A", "3.3.A")
  years <- c("y2000", "y2005")
  main <- c("crop", "past", "forestry", "primforest", "secdforest", "urban", "other")
  .new <- function(names, sets = NULL) new.magpie("A.1", years, names, fill = 0, sets = sets)

  startHr <- new.magpie(cells, "y1995", c("crop", "past", "range", "forestry", "primforest",
                                          "secdforest", "urban", "primother", "secdother"),
                        fill = c(1, 0.2, 0.5, 0.5, 1, 0.5, 0, 0.5, 0, 0.2, 0, 0.3,
                                 1, 0.5, 0, 0.5, 1, 1, 0.1, 0, 0.2, 0.5, 0.3, 1, 0.2, 0.5, 0.5))

  # cluster land (landAccounting), with the transitions that the land module fixes to 0
  dataLr <- list(land = .new(main),
                 forestry = .new(c("aff", "plant"), c("i", "j", "t", "type32")),
                 other = .new(c("othernat.ac0", "youngsecdf.ac0"), c("i", "j", "t", "othertype35", "ac")),
                 transitionsUpper = .new(paste(rep(main, each = 7), rep(main, 7), sep = ".")))
  dataLr[["transitionsUpper"]][, , ] <- Inf
  dataLr[["transitionsUpper"]][, , c(paste(setdiff(main, "primforest"), "primforest", sep = "."),
                                     "primforest.forestry", "primforest.other", "secdforest.other")] <- 0
  dataLr[["land"]][, "y2000", ] <- c(2.1, 2.2, 0.55, 1.4, 2.65, 0.3, 2.8)
  dataLr[["land"]][, "y2005", ] <- c(2.0, 2.3, 0.55, 1.4, 2.65, 0.3, 2.8)
  dataLr[["forestry"]][, , ] <- c(0.15, 0.15, 0.4, 0.4)
  dataLr[["other"]][, , ] <- c(2.57, 2.53, 0.23, 0.27)

  # the transitions by stage and sub-type (magpie4::landTransitions)
  stages <- c("disturbance.primforest.secdforest", "recoveryForestry.othernat.youngsecdf",
              "recovery.othernat.youngsecdf", "maturation.youngsecdf.secdforest", "optimisation.past.crop",
              "optimisation.othernat.crop", "optimisation.crop.plant", "optimisation.aff.othernat",
              "optimisation.youngsecdf.othernat", "optimisation.crop.past")
  transitionsLr <- new.magpie("A.1", years, stages, fill = 0, sets = c("i", "j", "t", "stage", "from", "to"))
  transitionsLr[, "y2000", stages[-c(2, 10)]] <- c(0.1, 0.1, 0.05, 0.3, 0.2, 0.1, 0.05, 0.02)
  transitionsLr[, "y2005", stages[c(2, 10)]] <- c(0.04, 0.1)

  cropSplitLr <- .new(c("crop_area", "crop_fallow", "crop_treecover"))
  cropSplitLr[, , ] <- c(6, 6, 2, 2, 1, 1) / 9 * 2
  cropAreaLr <- .new(c("tece.rainfed", "begr.irrigated"))
  cropAreaLr[, , ] <- c(0.75, 0.75, 0.25, 0.25)

  return(list(dataLr = dataLr, transitionsLr = transitionsLr, startHr = startHr,
              mapping = data.frame(cell = cells, cluster = "A.1"), cropSplitLr = cropSplitLr, cropAreaLr = cropAreaLr))
}

# The case with a second cluster B.2, a copy of cluster A.1 on the grid cells 4.4.B, 5.5.B and 6.6.B.
twoClusters <- function(case) {
  .copy <- function(x, items) {
    y <- x
    getItems(y, dim = 1, raw = TRUE) <- items
    return(mbind(x, y))
  }
  for (name in names(case[["dataLr"]])) {
    case[["dataLr"]][[name]] <- .copy(case[["dataLr"]][[name]], "B.2")
  }
  case[["transitionsLr"]] <- .copy(case[["transitionsLr"]], "B.2")
  case[["cropSplitLr"]] <- .copy(case[["cropSplitLr"]], "B.2")
  case[["cropAreaLr"]] <- .copy(case[["cropAreaLr"]], "B.2")
  case[["startHr"]] <- .copy(case[["startHr"]], c("4.4.B", "5.5.B", "6.6.B"))
  case[["mapping"]] <- rbind(case[["mapping"]], data.frame(cell = c("4.4.B", "5.5.B", "6.6.B"), cluster = "B.2"))
  return(case)
}

# The case with a transition of one stage ("stage.from.to") in a year, taken from the optimisation where the
# optimisation has that transition.
addStage <- function(case, item, year, value) {
  transitions <- case[["transitionsLr"]]
  optimisation <- sub("^[^.]*", "optimisation", item)
  if (optimisation %in% getNames(transitions)) {
    transitions[, year, optimisation] <- transitions[, year, optimisation] - value
  }
  added <- transitions[, , 1] * 0
  getNames(added) <- item
  added[, year, ] <- value
  case[["transitionsLr"]] <- mbind(transitions, added)
  return(case)
}

# Cluster sum of a cell transition in a year.
clusterFlow <- function(transitionsHr, year, item) sum(transitionsHr[, year, item])

# The sources of the calc functions for a case: config, clustermap, gdx symbols and 0.5 degree inputs, with
# available cropland and potential forest that do not bind unless given.
magpieSources <- function(case, config = list(), availableCropland = 100, potentialForest = 100) {
  cells <- case[["mapping"]][["cell"]]
  gms <- list(urban = "static", c29_marginal_land = "all", c35_pot_forest_correction = "on",
              c35_pot_forest_scenario = "cc", peatland = "off")
  gms[names(config)] <- config
  return(list(
    config = list(gms = gms),
    clustermap = case[["mapping"]],
    landAccounting = case[["dataLr"]],
    landTransitions = case[["transitionsLr"]],
    land = case[["cropSplitLr"]],
    crop = case[["cropAreaLr"]],
    landuseInit = case[["startHr"]],
    availableCropland = new.magpie(cells, NULL, "all", fill = availableCropland),
    potentialForest = new.magpie(cells, c("y2000", "y2005"), NULL, fill = potentialForest),
    snv = list(shareIso = c(A = 0, B = 0), fader = new.magpie("GLO", c("y2000", "y2005"), NULL, fill = 0),
               relocationX1 = 0.2, relocationX2 = 0.5),
    available = c("config", "landuseInit", "availableCropland", "potentialForest")
  ))
}

# Age classes of one cluster A with two grid cells and two years, age classes ac0, ac5, ac10 and acx: the
# model deforests 0.5 Mha of the oldest secondary forest in y2005 (grid cell c1 converts 0.5 Mha secondary
# forest to cropland), and grid cell c2 gains 0.2 Mha othernat.
ageCase <- function() {
  ac <- c("ac0", "ac5", "ac10", "acx")
  years <- c("y2000", "y2005")
  pools <- c("secdforest.total", "forestry.plant", "other.othernat", "other.youngsecdf")
  stocks <- new.magpie("A", years, paste(rep(pools, each = 4), ac, sep = "."), fill = 0,
                       sets = c("region", "year", "land", "type", "ac"))
  stocks[, "y2000", "secdforest.total"] <- c(0, 1, 1, 0)
  stocks[, "y2005", "secdforest.total"] <- c(0, 0, 1, 0.5)
  stocks[, "y2005", "other.othernat"] <- c(0.2, 0, 0, 0)

  processes <- c("disturbance.secdforest.total", "recovery.other.othernat", "maturation.other.youngsecdf",
                 paste(rep(c("outflow", "regrowth"), 4), rep(pools, each = 2), sep = "."))
  flows <- new.magpie("A", years, paste(rep(processes, each = 4), ac, sep = "."), fill = 0,
                      sets = c("region", "year", "process", "land", "type", "ac"))
  flows[, "y2005", "outflow.secdforest.total.acx"] <- 0.5

  cells <- c("c1", "c2")
  land <- new.magpie(cells, years, c("crop.total", pools), fill = 0)
  land[, , "secdforest.total"] <- c(1, 1, 0.5, 1)
  land[, "y2005", "other.othernat"] <- c(0, 0.2)
  transitions <- new.magpie(cells, years, c("optimisation.secdforest.total.crop.total",
                                            "optimisation.crop.total.other.othernat"), fill = 0)
  transitions["c1", "y2005", "optimisation.secdforest.total.crop.total"] <- 0.5
  transitions["c2", "y2005", "optimisation.crop.total.other.othernat"] <- 0.2

  return(list(sources = list(ageClasses = stocks, ageClassTransitions = flows,
                             clustermap = data.frame(cell = cells, cluster = "A")),
              land = list(landHr = land, transitionsHr = transitions)))
}
