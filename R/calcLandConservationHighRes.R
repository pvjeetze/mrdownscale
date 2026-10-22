#' calcLandConservationHighRes
#'
#' Land conservation of an input dataset by grid cell and land type: gridded conservation data fitted to
#' the conservation of the units of the input, by groups of land types.
#'
#' @details
#' input = "magpie": the data are selected by the switches of module 22: base protection
#' (\code{c22_base_protect}, \code{c22_base_protect_noselect}) by year as \code{p22_wdpa_baseline}, the
#' last year of the data after it, for the policy countries and the others, plus future conservation
#' in the conservation priority areas (\code{c22_protect_scenario}, \code{c22_protect_scenario_noselect})
#' with the fader of the run, for natural vegetation only. The conservation of the model is
#' \code{pm_land_conservation} (protection and restoration), fitted in the groups pasture and natural
#' vegetation (primary forest, secondary forest and other land, as a whole because the land types of
#' the data need not match those of the start map).
#'
#' Fit, for each group: without land (\code{onLand = FALSE}), the data are scaled with one factor per
#' unit to the unit conservation of the input. On the land (\code{onLand = TRUE}, the land of
#' \code{\link{calcLandTransitionsHighRes}}), the unit conservation goes first to the grid cells with
#' conservation in the data, in proportion to it, up to their land, then to the other grid cells in
#' proportion to their land (\code{\link{toolAllocateByFilling}}); it is split over the land types of the
#' group as the land of the grid cell. Land types outside the groups keep the data without land and are 0
#' on the land.
#'
#' @param input name of an input dataset, currently only "magpie"
#' @param onLand FALSE to scale the data, TRUE to fit the conservation on the downscaled land
#' @param harmonization harmonisation of the downscaled land with \code{onLand = TRUE}, see
#' \code{\link{calcLandTransitionsHighRes}}
#' @return magclass object of land conservation by grid cell and land type (Mha)
#' @author Patrick v. Jeetze
#' @seealso \code{\link{calcLandConstraints}}
calcLandConservationHighRes <- function(input, onLand = FALSE, harmonization = "none") {
  # Sum of the types of each group (grid cells or units x years x groups).
  .groupSums <- function(x, groups) {
    out <- list()
    for (group in names(groups)) {
      out[[group]] <- setNames(dimSums(x[, , groups[[group]]], dim = 3), group)
    }
    return(mbind(out))
  }

  # =================================================================================================
  # Conservation data on the grid and conservation of the units, by input
  # =================================================================================================

  if (input == "magpie") {
    # data selected by the switches of module 22
    cfg <- readSource("MagpieRunFiles", subtype = "config")
    mapping <- readSource("MagpieFulldataGdx", subtype = "clustermap")
    conservation <- readSource("MagpieFulldataGdx", subtype = "landConservation")
    wdpaHr <- readSource("MagpieRunFiles", subtype = "wdpa")
    landConservationLr <- conservation[["landConservation"]]
    cells <- as.character(mapping[["cell"]])
    years <- getYears(landConservationLr)

    dataHr <- new.magpie(cells, years, getItems(wdpaHr, dim = 3.2), fill = 0, sets = c("x.y.iso", "year", "data"))
    consvIso <- conservation[["policyCountries"]]
    consvIso <- consvIso[consvIso %in% getItems(wdpaHr, dim = 1.3)]
    if (length(consvIso) == 0) {
      warning("No countries selected in the land conservation downscaling. Results may be erroneous")
    }

    # base protection by year as p22_wdpa_baseline, the last year of the data after it
    baseProtect <- cfg[["gms"]][["c22_base_protect"]]
    baseProtectNoselect <- cfg[["gms"]][["c22_base_protect_noselect"]]
    if (baseProtectNoselect != "none") {
      dataHr[, , ] <- collapseDim(wdpaHr[, nyears(wdpaHr), baseProtectNoselect], dim = 3.1)
      dataHr[, getYears(wdpaHr), ] <- collapseDim(wdpaHr[, , baseProtectNoselect], dim = 3.1)
    }
    if (baseProtect != "none") {
      dataHr[consvIso, , ] <- collapseDim(wdpaHr[consvIso, nyears(wdpaHr), baseProtect], dim = 3.1)
      dataHr[consvIso, getYears(wdpaHr), ] <- collapseDim(wdpaHr[consvIso, , baseProtect], dim = 3.1)
    } else {
      dataHr[consvIso, , ] <- 0
    }

    # future conservation in the conservation priority areas, with the fader
    protectScenario <- cfg[["gms"]][["c22_protect_scenario"]]
    protectScenarioNoselect <- cfg[["gms"]][["c22_protect_scenario_noselect"]]
    if (!all(c(protectScenario, protectScenarioNoselect) %in% "none")) {
      if ("conservationPriority" %in% readSource("MagpieRunFiles", subtype = "available")) {
        consvPrioAll <- readSource("MagpieRunFiles", subtype = "conservationPriority")
        consvPrioHr <- new.magpie(cells_and_regions = cells, names = getNames(consvPrioAll, dim = 2),
                                  fill = 0, sets = c("x.y.iso", "year", "data"))
        if (protectScenarioNoselect != "none") {
          consvPrioHr <- collapseDim(consvPrioAll[, , protectScenarioNoselect], dim = 3.1)
        }
        if (protectScenario != "none") {
          consvPrioHr[consvIso, , ] <- collapseDim(consvPrioAll[consvIso, , protectScenario], dim = 3.1)
        } else {
          consvPrioHr[consvIso, , ] <- 0
        }
        # future conservation only pertains to natural vegetation
        consvPrioHr[, , c("crop", "past", "forestry", "urban")] <- 0
        dataHr <- dataHr + consvPrioHr * conservation[["fader"]][, years, ]
      } else {
        warning("Future land conservation used in MAgPIE run but high resolution ",
                "conservation priority data for the downscaling not found.")
      }
    }

    # land conservation of the model (protection and restoration), fitted in the groups pasture and
    # natural vegetation
    modelLr <- dimSums(landConservationLr, dim = 3.2)
    groups <- list(past = "past", natveg = c("primforest", "secdforest", "other"))
  } else {
    stop("Unsupported input dataset \"", input, "\"")
  }

  # =================================================================================================
  # Fit to the conservation of the units
  # =================================================================================================

  modelGroupLr <- .groupSums(modelLr, groups)
  dataGroupHr <- .groupSums(dataHr, groups)
  conservationHr <- dataHr
  if (!onLand) {
    # the data scaled with one factor per unit
    scalingHr <- modelGroupLr / toolAggregate(dataGroupHr, mapping, from = "cell", to = "cluster")
    scalingHr[is.na(scalingHr) | is.infinite(scalingHr)] <- 1
    scalingHr <- toolAggregate(scalingHr, mapping, from = "cluster", to = "cell")
    for (group in names(groups)) {
      conservationHr[, , groups[[group]]] <- scalingHr[, , group] * dataHr[, , groups[[group]]]
    }
  } else {
    # land conservation achieved on the downscaled land
    landHr <- calcOutput("LandTransitionsHighRes", input = input, harmonization = harmonization,
                         aggregate = FALSE)[["landHr"]]
    landHr <- dimSums(landHr, dim = 3.2)
    stopifnot("the downscaled land lacks grid cells of the mapping" = all(cells %in% getCells(landHr)))
    landHr <- landHr[cells, years, ]
    landGroupHr <- .groupSums(landHr, groups)
    getSets(landGroupHr) <- getSets(dataGroupHr)

    # the conservation of each unit: first in the grid cells with conservation in the data, in proportion to
    # it, up to their land, then the rest in the other grid cells, in proportion to their land
    protectedGroupHr <- toolAllocateByFilling(modelGroupLr, pmin(dataGroupHr, landGroupHr), landGroupHr, mapping)

    # the conservation of a group split over its types as the land of the grid cell
    conservationHr[, , ] <- 0
    for (group in names(groups)) {
      types <- groups[[group]]
      share <- landHr[, , types] / dimSums(landHr[, , types], dim = 3)
      share[is.na(share)] <- 0
      conservationHr[, , types] <- collapseDim(protectedGroupHr[, , group], dim = 3) * share
    }
  }

  return(list(x = conservationHr, isocountries = FALSE, unit = "Mha", min = 0,
              description = "Land conservation of the input by grid cell and land type"))
}
