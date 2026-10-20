#' calcLandHighRes
#'
#' This function performs the downscaling: It calculates a high resolution dataset
#' from the low resolution input dataset and the high resolution target dataset
#' using the given downscaling method.
#'
#' @details With \code{downscaling = "transitions"} and \code{harmonization = "none"}, the land is
#' downscaled along the transitions of the input from its own start map
#' (\code{\link{calcLandTransitionsHighRes}}), without harmonisation and recategorisation: the land types
#' of the input are the categories of the target (only \code{target = "landuseinit"}), the grid cells
#' those of the input (names "x.y.iso"), and the start is labelled y1985. Otherwise the harmonised land
#' is downscaled from the target map at the end of history.
#'
#' @inheritParams calcLandInput
#' @param target name of a target dataset
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @param yearsSubset vector of years to keep in the output dataset; NULL keeps all years (only with
#' \code{downscaling = "transitions"} and \code{harmonization = "none"})
#' @param harmonization name of harmonization method, see \code{\link{toolGetHarmonizer}}, or "none"
#' (only with \code{downscaling = "transitions"})
#' @param downscaling name of downscaling method: "proportional" (\code{luscale::interpolate2}) or "transitions"
#' (\code{\link{calcLandTransitionsHighRes}}, or \code{\link{toolDownscaleTransitions}} on harmonised states)
#' @return downscaled land use data
#' @author Jan Philipp Dietrich, Pascal Sauer, Patrick v. Jeetze
calcLandHighRes <- function(input, target, harmonizationPeriod, yearsSubset, harmonization, downscaling) {
  if (downscaling == "transitions" && harmonization == "none") {
    if (target != "landuseinit") {
      stop("downscaling = \"transitions\" with harmonization = \"none\" supports only target = \"landuseinit\"")
    }
    landHr <- calcOutput("LandTransitionsHighRes", input = input, target = target, harmonization = harmonization,
                         aggregate = FALSE)[["landHr"]]
    out <- dimSums(landHr, dim = 3.2)
    if (!is.null(yearsSubset)) {
      out <- out[, getYears(out, as.integer = TRUE) %in% c(1985, yearsSubset), ]
    }

    targetCategories <- unique(toolGetMapping("referenceMappings/landuseinit.csv", where = "mrdownscale")[["data"]])
    toolExpectTrue(setequal(getItems(out, dim = 3), targetCategories), "Land categories match target definition")
    toolExpectTrue(all(out >= 0), "All values are >= 0")
    outSum <- dimSums(out, dim = 3)
    toolExpectLessDiff(outSum, outSum[, 1, ], 10^-5, "Total land area per cell in output stays constant over time")
    return(list(x = out, class = "magpie", isocountries = FALSE, unit = "Mha", min = 0,
                description = "Land use data downscaled along the transitions of the input"))
  }

  x <- calcOutput("LandHarmonized", input = input, target = target,
                  harmonizationPeriod = harmonizationPeriod, harmonization = harmonization,
                  aggregate = FALSE)
  x <- x[, getYears(x, as.integer = TRUE) %in% yearsSubset, ]

  hp1 <- harmonizationPeriod[1]

  xTarget <- calcOutput("LandTarget", target = target, endOfHistory = hp1, aggregate = FALSE)
  xTarget <- as.magpie(xTarget[[terra::time(xTarget) %in% yearsSubset]])

  histYears <- intersect(getYears(x, as.integer = TRUE),
                         getYears(xTarget, as.integer = TRUE))

  x <- x[, , getItems(xTarget, 3)]

  landTargetLowRes <- calcOutput("LandTargetLowRes", input = input, target = target,
                                 endOfHistory = hp1, aggregate = FALSE)
  landTargetLowRes <- landTargetLowRes[, max(histYears), getItems(x, 3)]
  stopifnot(setequal(getItems(landTargetLowRes, 1), getItems(x, 1)))

  mapping <- calcOutput("ResolutionMapping", input = input, target = target, aggregate = FALSE)

  if (downscaling == "proportional" || downscaling == "magpieClassic") {
    out <- toolDownscaleSimpleProportional(x[, getYears(x, as.integer = TRUE) >= max(histYears), ],
                                           xTarget[, max(histYears), ],
                                           xTargetLowRes = landTargetLowRes,
                                           mapping = mapping)
    out <- out[, getYears(out, as.integer = TRUE) > max(histYears), ]
  } else if (downscaling == "transitions") {
    out <- toolDownscaleTransitions(x[, getYears(x, as.integer = TRUE) >= max(histYears), ],
                                    xTarget[, max(histYears), ],
                                    mapping = mapping)
  } else {
    stop("Unsupported downscaling method \"", downscaling, "\"")
  }

  out <- mbind(xTarget[, histYears, ], out)

  out <- toolReplaceExpansion(out, "primf", "secdf", warnThreshold = 100)
  out <- toolReplaceExpansion(out, "primn", "secdn", warnThreshold = 100)

  toolExpectTrue(identical(unname(getSets(out)), c("x", "y", "year", "data")),
                 "Dimensions are named correctly")
  toolExpectTrue(setequal(getItems(out, dim = 3), getItems(x, dim = 3)),
                 "Land categories remain unchanged")
  toolExpectLessDiff(out[, max(histYears), ], xTarget[, max(histYears), ], 10^-5,
                     paste("In", max(histYears), "output equals target"))
  toolExpectTrue(all(out >= 0), "All values are >= 0")

  outSum <- dimSums(out, dim = 3)
  toolExpectLessDiff(outSum, outSum[, 1, ], 10^-5,
                     "Total land area per cell in output stays constant over time")

  globalSumIn <- dimSums(x[, getYears(out), ], dim = 1)
  globalSumOut <- dimSums(out, dim = 1)
  toolExpectLessDiff(dimSums(globalSumIn, 3), dimSums(globalSumOut, 3), 10^-5,
                     "Total global land area remains unchanged")
  toolExpectLessDiff(globalSumIn, globalSumOut, 10^-5,
                     "Global area of each land type remains unchanged")
  toolPrimExpansionCheck(out)

  return(list(x = out,
              class = "magpie",
              isocountries = FALSE,
              unit = "Mha",
              min = 0,
              description = "Downscaled land use data"))
}
