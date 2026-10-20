#' calcLandConstraints
#'
#' Constraints on the grid cell land of an input dataset for the downscaling along its transitions
#' (\code{\link{calcLandTransitionsHighRes}}): bounds on the summed area of groups of land types per grid
#' cell, each with a priority. Where they cannot all hold within a unit, the one with the lowest
#' priority is exceeded first (see \code{\link{toolAllocateTransitions}}).
#'
#' @details
#' input = "magpie": the constraints of the model, which apply in all stages:
#' \itemize{
#'   \item potentialForest: forest establishment within the potential forest area
#'   (\code{q35_max_forest_establishment}), potential forest by \code{c35_pot_forest_correction} and
#'   \code{c35_pot_forest_scenario} as in modules/35_natveg/pot_forest_may24/input.gms;
#'   \item crop: cropland within the available cropland (\code{c29_marginal_land}) reduced by the SNV
#'   share and its fader (\code{q29_avl_cropland});
#'   \item natveg and past: protected natural vegetation (as a whole, because the land types of the
#'   protected-area data need not match those of the start map) and pasture (\code{vm_land.lo}) from
#'   \code{\link{calcLandConservationHighRes}}, if the run has land conservation data.
#' }
#' Priorities: protection 1, potentialForest 2, crop 3 (protection is exceeded first).
#'
#' @param input name of an input dataset, currently only "magpie"
#' @return Named list of constraints, each a list with \code{types} (land types; main types stand for
#' their sub-types), \code{upper} or \code{lower} (magclass object, grid cells x years) and
#' \code{priority}
#' @author Patrick v. Jeetze
#' @seealso \code{\link{calcLandTransitionsInput}}
calcLandConstraints <- function(input) {
  if (input == "magpie") {
    priorities <- c(protection = 1, potentialForest = 2, crop = 3)
    cfg <- readSource("MagpieRunFiles", subtype = "config")

    # potential forest area of the run: corrected for grassland ecoregions or not, and its climate
    # change scenario
    if (identical(cfg[["gms"]][["c35_pot_forest_correction"]], "off")) {
      potentialForest <- readSource("MagpieRunFiles", subtype = "potentialForestUncorrected")
    } else {
      potentialForest <- readSource("MagpieRunFiles", subtype = "potentialForest")
    }
    years <- getYears(potentialForest, as.integer = TRUE)
    sourceYears <- years
    if (identical(cfg[["gms"]][["c35_pot_forest_scenario"]], "nocc")) {
      sourceYears[] <- 1995
    } else if (identical(cfg[["gms"]][["c35_pot_forest_scenario"]], "nocc_hist")) {
      fixYear <- readSource("MagpieFulldataGdx", subtype = "climateFixYear")
      sourceYears[years > fixYear] <- fixYear
    }
    potentialForest <- setYears(potentialForest[, paste0("y", sourceYears), ], getYears(potentialForest))

    # available cropland reduced by the SNV share of the country of each grid cell and its fader
    croplandHr <- readSource("MagpieRunFiles", subtype = "availableCropland")[, , cfg[["gms"]][["c29_marginal_land"]]]
    snv <- readSource("MagpieFulldataGdx", subtype = "snv")
    snvShareHr <- new.magpie(getCells(croplandHr),
                             fill = unname(snv[["shareIso"]][getItems(croplandHr, dim = 1.3, full = TRUE)]))
    snvFader <- collapseDim(snv[["fader"]], dim = 1)

    constraints <- list(
      potentialForest = list(types = c("forestry", "primforest", "secdforest", "youngsecdf"),
                             upper = potentialForest, priority = priorities[["potentialForest"]]),
      crop = list(types = "crop", upper = croplandHr * (1 - snvShareHr * snvFader), priority = priorities[["crop"]])
    )

    # protected natural vegetation and pasture
    if ("wdpa" %in% readSource("MagpieRunFiles", subtype = "available")) {
      protectionHr <- calcOutput("LandConservationHighRes", input = input, onLand = FALSE, aggregate = FALSE)
      natveg <- c("primforest", "secdforest", "other")
      constraints[["natveg"]] <- list(types = natveg, priority = priorities[["protection"]],
                                      lower = dimSums(protectionHr[, , natveg], dim = 3))
      constraints[["past"]] <- list(types = "past", lower = protectionHr[, , "past"],
                                    priority = priorities[["protection"]])
    }
  } else {
    stop("Unsupported input dataset \"", input, "\"")
  }

  return(list(x = constraints, class = "list", unit = "Mha",
              description = "Constraints on the grid cell land with their priorities"))
}
