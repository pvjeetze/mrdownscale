#' calcLandReportBtC2
#'
#' The land of the downscaling along the transitions (\code{\link{calcLandTransitionsHighRes}}) in the land
#' classes of the project Bending the Curve 2.0 (BtC2), by grid cell. The classes only reclassify the
#' downscaled land (\code{cell.land_0.5.mz}); the land itself is not changed.
#'
#' @details
#' Classes of the state of the land, from the land split (\code{\link{calcLandReport}} with
#' \code{outputFormat = "downscaledmagpieSplit"}):
#' \itemize{
#'   \item cropland_other and cropland_2Gbioen: cropland without and with the crop area of bioenergy crops;
#'   \item grassland: pasture and rangeland;
#'   \item nat_regen_forest: primary and secondary forest, NDC forest, and afforestation if it is not
#'   established as plantation (\code{s32_aff_plantation});
#'   \item forest_planted: timber plantations, other planted forest, and afforestation as plantation;
#'   \item built_up_areas: urban land.
#' }
#' Classes of the history of other land (othernat and youngsecdf), along the grid cell transitions from the
#' first year on (\code{\link{toolGriddedLandOrigin}}): area that came from cropland (rest_abn_crop), pasture
#' (rest_abn_grass), or timber plantations and other planted forest (rest_forest) keeps its class until it
#' leaves other land; the rest (other land of the first year, from NDC forest, afforestation and urban land)
#' is other. Maturation to secondary forest and the conversion to NDC forest and to afforestation that is
#' not a plantation take the area that entered other land first; all other conversions the area that
#' entered last.
#'
#' Restoration of the optimisation (the lower bounds \code{p35_land_restoration}): in each cluster the bound
#' as share of the transitions into the restored type in the optimisation stage; the transitions of each
#' grid cell are split by this share. All outflows of other land except recovery and maturation (conversions
#' and harvest) take restored land only after all other area of the type in the grid cell; where the land
#' loses more, restored land is taken. Where the same happens on the cluster sums, the model converts the
#' restored land itself (MAgPIE protects restored land only as area of its type in the cluster), reported as a
#' note; restored land taken beyond that comes from the placement in grid cells and is reported with a warning.
#'
#' Land that enters secondary forest directly (from cropland, pasture, forestry or urban land) counts as
#' recovering land of its origin until it reaches the age at which MAgPIE matures youngsecdf: the first age
#' class with a vegetation carbon density of secondary forest above 20 tC/ha in its cluster
#' (modules/35_natveg/pot_forest_may24/presolve.gms). A loss of secondary forest in a grid cell takes from
#' this land in proportion to the area of secondary forest, from restored land only beyond all other
#' secondary forest.
#'
#' @param input name of an input dataset, currently only "magpie"
#' @return land by grid cell and class (Mha) in the years of the input
#' @author Patrick v. Jeetze
#' @seealso \code{\link{fullBTC2}}, \code{\link{toolGriddedLandOrigin}}
calcLandReportBtC2 <- function(input = "magpie") {
  # Land report of the run in one of the MAgPIE formats, with the arguments of fullDOWNSCALEDMAGPIE (same cache).
  .report <- function(outputFormat) {
    x <- calcOutput("LandReport", outputFormat = outputFormat, input = input, harmonizationPeriod = c(2015, 2050),
                    yearsSubset = seq(1995, 2100, 5), harmonization = "none", downscaling = "transitions",
                    aggregate = FALSE)
    return(x[, "y1985", , invert = TRUE])
  }

  # Splits the transitions into secondary forest and other land into the restored types by the restoration
  # share of their cluster (restoration bound per transitions of the cluster); the restored part has the source
  # "<source>Restored". The bounds apply to all transitions of the model, which magpie4::landTransitions
  # splits into the stages snvRelocation, peatRewetting and optimisation.
  .splitRestoration <- function(transitionsHr, restorationLr, mapping) {
    stage <- getItems(transitionsHr, "stage", full = TRUE)
    fromMain <- getItems(transitionsHr, "from", full = TRUE)
    toMain <- getItems(transitionsHr, "to", full = TRUE)
    sources <- list(secdforest = c("crop", "past", "forestry"), other = c("crop", "past"))
    years <- getYears(transitionsHr)
    restored <- list()
    for (k in names(sources)) {
      select <- which(stage %in% c("snvRelocation", "peatRewetting", "optimisation") & toMain == k &
                        fromMain %in% sources[[k]])
      flowLr <- toolAggregate(dimSums(transitionsHr[, , select], dim = 3), mapping, from = "cell", to = "cluster")
      shareLr <- pmin(restorationLr[getItems(flowLr, dim = 1), years, k] / setNames(flowLr, k), 1)
      shareLr[is.na(shareLr) | is.infinite(shareLr)] <- 0
      shareHr <- toolAggregate(shareLr, mapping, from = "cluster", to = "cell")[getCells(transitionsHr), , ]
      restored[[k]] <- transitionsHr[, , select] * setNames(shareHr, NULL)
      getNames(restored[[k]]) <- paste(stage[select], paste0(fromMain[select], "Restored"),
                                       getItems(transitionsHr, "fromType", full = TRUE)[select], toMain[select],
                                       getItems(transitionsHr, "toType", full = TRUE)[select], sep = ".")
      transitionsHr[, , select] <- transitionsHr[, , select] * setNames(1 - shareHr, NULL)
    }
    return(mbind(transitionsHr, restored[["secdforest"]], restored[["other"]]))
  }

  # Land that entered secondary forest directly, by grid cell and origin (sources in restoredFrom), until its
  # maturation age (grid cells x years of transitionsHr); secondary forest lost by a grid cell takes from it
  # in proportion, from the protected origins only beyond all other secondary forest. The cohorts are an array:
  # slices are written in a loop, which in magclass takes about 4000 times as long (see toolGriddedLandOrigin).
  .directSecdforest <- function(transitionsHr, secdforestStartHr, restoredFrom, protected, maturityAge) {
    stage <- getItems(transitionsHr, "stage", full = TRUE)
    from <- getItems(transitionsHr, "from", full = TRUE)
    to <- getItems(transitionsHr, "to", full = TRUE)
    years <- getYears(transitionsHr, as.integer = TRUE)
    nCell <- ncells(transitionsHr)
    secdforest <- as.vector(secdforestStartHr)
    cohorts <- array(0, c(nCell, length(restoredFrom), length(years)))
    out <- new.magpie(getCells(transitionsHr), years, restoredFrom, fill = 0)
    lost <- stats::setNames(rep(0, length(years)), years)
    for (i in seq_along(years)) {
      for (e in seq_len(i - 1)) {
        cohorts[years[i] - years[e] >= maturityAge[, i], , e] <- 0
      }
      flows <- matrix(transitionsHr[, i, ], nrow = nCell)
      for (s in unique(stage)) {
        outflow <- rowSums(flows[, stage == s & from == "secdforest" & to != "secdforest", drop = FALSE])
        inflow <- rowSums(flows[, stage == s & to == "secdforest" & from != "secdforest", drop = FALSE])
        protectedArea <- rowSums(cohorts[, protected, , drop = FALSE])
        free <- secdforest - protectedArea
        cohorts[, !protected, ] <- cohorts[, !protected, , drop = FALSE] *
          (1 - pmin(outflow / replace(free, free <= 0, 1), 1))
        excess <- pmin(pmax(outflow - free, 0), protectedArea)
        cohorts[, protected, ] <- cohorts[, protected, , drop = FALSE] *
          (1 - excess / replace(protectedArea, protectedArea <= 0, 1))
        lost[i] <- lost[i] + sum(excess)
        # the outflow beyond the area at the start of the stage passes through from the inflow of the stage
        remaining <- 1 - pmin(pmax(outflow - secdforest, 0) / replace(inflow, inflow <= 0, 1), 1)
        for (k in seq_along(restoredFrom)) {
          entering <- stage == s & to == "secdforest" & from == restoredFrom[k]
          cohorts[, k, i] <- cohorts[, k, i] + rowSums(flows[, entering, drop = FALSE]) * remaining
        }
        secdforest <- secdforest + inflow - outflow
      }
      out[, i, ] <- rowSums(cohorts, dims = 2)
    }
    attr(out, "lost") <- lost
    return(out)
  }

  if (input != "magpie") {
    stop("Can only report for input = magpie")
  }
  landHr <- .report("downscaledmagpie")
  splitHr <- .report("downscaledmagpieSplit")
  transitionsHr <- calcOutput("LandTransitionsHighRes", input = input, aggregate = FALSE)[["transitionsHr"]]
  transitionsHr <- transitionsHr[, getYears(landHr)[-1], ]
  mapping <- readSource("MagpieFulldataGdx", subtype = "clustermap")
  cfg <- readSource("MagpieRunFiles", subtype = "config")

  # =================================================================================================
  # Classes of the state of the land
  # =================================================================================================

  splitItems <- getItems(splitHr, dim = 3)
  classOf <- stats::setNames(rep(NA_character_, length(splitItems)), splitItems)
  classOf[startsWith(splitItems, "crop.")] <- "cropland_other"
  classOf[startsWith(splitItems, "crop.kbe.")] <- "cropland_2Gbioen"
  classOf[startsWith(splitItems, "past.")] <- "grassland"
  classOf[startsWith(splitItems, "primforest.") | startsWith(splitItems, "secdforest.")] <- "nat_regen_forest"
  classOf[startsWith(splitItems, "forestry.")] <- "forest_planted"
  classOf["forestry.ndc.total"] <- "nat_regen_forest"
  affPlantation <- as.logical(cfg$gms$s32_aff_plantation)
  classOf["forestry.aff.total"] <- if (affPlantation) "forest_planted" else "nat_regen_forest"
  classOf[startsWith(splitItems, "urban.")] <- "built_up_areas"
  stateClasses <- c("cropland_other", "cropland_2Gbioen", "grassland", "nat_regen_forest", "forest_planted",
                    "built_up_areas")
  stateHr <- mbind(lapply(stateClasses, function(class) {
    setNames(dimSums(splitHr[, , names(which(classOf == class))], dim = 3), class)
  }))

  # =================================================================================================
  # Classes of the history of other land
  # =================================================================================================

  # restored land as source "<source>Restored"; labels "stage.from.to" with the sub-types of other land and
  # forestry as destination, and as source managed forest (timber plantations, other planted forest) as
  # "timber", NDC forest and afforestation as "forestry"
  transitionsHr <- .splitRestoration(transitionsHr, readSource("MagpieFulldataGdx", subtype = "restoration"),
                                     mapping)
  # cropland, pasture and forestry become secondary forest only by restoration (q35_secdforest_regeneration), so
  # none of their transitions into secondary forest stays without the restored source
  unrestored <- getItems(transitionsHr, "to", full = TRUE) == "secdforest" &
    getItems(transitionsHr, "from", full = TRUE) %in% c("crop", "past", "forestry")
  toolExpectTrue(sum(transitionsHr[, , unrestored]) <= 1e-5,
                 "Transitions from cropland, pasture and forestry into secondary forest are restoration")
  restored <- endsWith(getItems(transitionsHr, "from", full = TRUE), "Restored")
  fromLabel <- getItems(transitionsHr, "fromType", full = TRUE)
  fromLabel[fromLabel == "total"] <- getItems(transitionsHr, "from", full = TRUE)[fromLabel == "total"]
  fromLabel[fromLabel %in% c("plant", "other_planted")] <- "timber"
  fromLabel[fromLabel %in% c("aff", "ndc")] <- "forestry"
  fromLabel[restored] <- paste0(sub("Restored$", "", fromLabel[restored]), "Restored")
  toLabel <- getItems(transitionsHr, "toType", full = TRUE)
  toLabel[toLabel == "total"] <- getItems(transitionsHr, "to", full = TRUE)[toLabel == "total"]
  label <- paste(getItems(transitionsHr, "stage", full = TRUE), fromLabel, toLabel, sep = ".")
  transitionsHr <- mbind(lapply(unique(label), function(l) {
    setNames(dimSums(transitionsHr[, , which(label == l)], dim = 3), l)
  }))
  getSets(transitionsHr, fulldim = FALSE)[3] <- "stage.from.to"

  restoredOther <- c("cropRestored", "pastRestored")
  naturalForestry <- c("ndc", if (!affPlantation) "aff")
  oldest <- c("youngsecdf.secdforest", paste(rep(c("othernat", "youngsecdf"), each = length(naturalForestry)),
                                             naturalForestry, sep = "."))
  outflowOrder <- c(default = "newest", stats::setNames(rep("oldest", length(oldest)), oldest))
  last <- list(default = restoredOther, othernat.youngsecdf = character(0), youngsecdf.secdforest = character(0))
  startHr <- setNames(splitHr[, 1, c("other.othernat.total", "other.youngsecdf.total")], c("othernat", "youngsecdf"))
  originHr <- toolGriddedLandOrigin(transitionsHr, startHr, c("othernat", "youngsecdf"), order = outflowOrder,
                                    last = last)

  # land that entered secondary forest directly, as recovering land until its maturation age
  densityLr <- readSource("MagpieFulldataGdx", subtype = "secdforestCarbonDensity")[, getYears(transitionsHr), ]
  # age of each age class from its name (ac0, ac5, ...; acx, the oldest, never matures)
  age <- suppressWarnings(as.numeric(sub("^ac", "", getItems(densityLr, dim = 3))))
  # maturation threshold of 20 tC/ha, keep in sync with p35_maturesecdf in
  # modules/35_natveg/pot_forest_may24/presolve.gms
  maturityAgeLr <- apply(as.array(densityLr > 20), c(1, 2), function(x) min(age[x & !is.na(age)], Inf))
  cellCluster <- mapping[["cluster"]][match(getCells(landHr), mapping[["cell"]])]
  restoredFrom <- c("crop", "past", "timber", "forestry", "urban", "cropRestored", "pastRestored",
                    "timberRestored", "forestryRestored")
  directHr <- .directSecdforest(transitionsHr, landHr[, 1, "secdforest"], restoredFrom,
                                endsWith(restoredFrom, "Restored"), maturityAgeLr[cellCluster, , drop = FALSE])
  takenRestored <- attr(originHr, "takenLast") + attr(directHr, "lost")
  directHr <- mbind(setYears(directHr[, 1, ] * 0, getYears(landHr)[1]), directHr)

  # the same on the cluster sums: restored land that the model converts itself (MAgPIE protects restored land
  # only as area of its type in the cluster); the rest of the restored land taken comes from the placement
  transitionsLr <- toolAggregate(transitionsHr, mapping, from = "cell", to = "cluster")
  originLr <- toolGriddedLandOrigin(transitionsLr, toolAggregate(startHr, mapping, from = "cell", to = "cluster"),
                                    c("othernat", "youngsecdf"), order = outflowOrder, last = last)
  directLr <- .directSecdforest(transitionsLr,
                                toolAggregate(landHr[, 1, "secdforest"], mapping, from = "cell", to = "cluster"),
                                restoredFrom, endsWith(restoredFrom, "Restored"),
                                maturityAgeLr[getCells(transitionsLr), , drop = FALSE])
  takenByModel <- attr(originLr, "takenLast") + attr(directLr, "lost")
  takenByPlacement <- pmax(takenRestored - takenByModel, 0)

  # the history classes, each of its origins in other land and in direct secondary forest
  originsOf <- list(rest_abn_crop = c("crop", "cropRestored"), rest_abn_grass = c("past", "pastRestored"),
                    rest_forest = c("timber", "timberRestored"))
  recoveringHr <- mbind(lapply(names(originsOf), function(class) {
    items <- intersect(outer(c("othernat", "youngsecdf"), originsOf[[class]], paste, sep = "."),
                       getItems(originHr, dim = 3))
    area <- dimSums(directHr[, , originsOf[[class]]], dim = 3)
    if (length(items) > 0) {
      area <- dimSums(originHr[, , items], dim = 3) + area
    }
    setNames(area, class)
  }))
  otherHr <- setNames(landHr[, , "other"] + dimSums(directHr, dim = 3) - dimSums(recoveringHr, dim = 3), "other")
  stateHr[, , "nat_regen_forest"] <- stateHr[, , "nat_regen_forest"] - dimSums(directHr, dim = 3)

  if (max(takenByModel) > 1e-6) {
    toolStatusMessage("note", paste0("MAgPIE converted restored land (it protects restored land only as area of ",
                                     "its type in the cluster), up to ", signif(max(takenByModel), 3),
                                     " Mha in a time step (", names(which.max(takenByModel)), ")"))
  }
  if (max(takenByPlacement) > 1e-6) {
    toolStatusMessage("warn", paste0("The grid cell placement converted restored land beyond the model, up to ",
                                     signif(max(takenByPlacement), 3), " Mha in a time step (",
                                     names(which.max(takenByPlacement)), ")"))
  }

  x <- mbind(stateHr[, , setdiff(stateClasses, "built_up_areas")], recoveringHr, otherHr,
             stateHr[, , "built_up_areas"])
  getSets(x, fulldim = FALSE) <- c("x.y.iso", "year", "class")
  toolExpectLessDiff(dimSums(x, dim = 3), dimSums(landHr, dim = 3), 1e-6,
                     "The classes cover all land of the grid cells")
  toolExpectTrue(min(x) >= -1e-6, "No class has negative area")
  return(list(x = x, isocountries = FALSE, unit = "Mha", cache = FALSE,
              description = "MAgPIE land by grid cell in the land classes of Bending the Curve 2.0"))
}
