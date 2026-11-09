#' calcLandReportScenarioMIP
#'
#' Convert the downscaled land use data to the format required by ScenarioMIP.
#'
#' @inheritParams calcLandInput
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @param yearsSubset vector of years to keep in the output dataset
#' @param harmonization name of harmonization method, see \code{\link{toolGetHarmonizer}}
#' @param downscaling name of downscaling method: "proportional" (\code{luscale::interpolate2}) or "transitions"
#' (\code{\link{toolDownscaleTransitions}})
#' @return land use data
#' @author Pascal Sauer
calcLandReportScenarioMIP <- function(input, harmonizationPeriod, yearsSubset, harmonization, downscaling) {
  landHighRes <- calcOutput("LandHighRes", input = input, target = "luh3",
                            harmonizationPeriod = harmonizationPeriod, yearsSubset = yearsSubset,
                            harmonization = harmonization, downscaling = downscaling,
                            aggregate = FALSE)
  cellArea <- readSource("LUH3", subtype = "cellArea", convert = FALSE)
  cellArea <- collapseDim(as.magpie(cellArea), 3)

  cropData <- toolCropData(landHighRes, cellArea)

  nonCropData <- landHighRes[, , c("primf", "primn", "secdf", "secdn", "urban", "pastr", "range", "pltns")]
  nonCropData <- nonCropData * 10000 / cellArea[getItems(cellArea, 1) %in% getItems(nonCropData, 1), , ]

  out <- mbind(cropData, nonCropData)

  if (max(out) > 1.0001) {
    toolStatusMessage("warn", paste0("Some shares are > 1 (max: ", max(out), "), setting those to 1"))
  }
  out[out > 1] <- 1

  toolExpectTrue(min(out) >= 0, "All values are >= 0")
  toolExpectTrue(max(out) <= 1, "All shares are <= 1")
  landShares <- c("c3ann", "c3nfx",  "c3per", "c4ann", "c4per", "primf",
                  "primn", "secdf", "secdn", "urban", "pastr", "range", "pltns")
  toolExpectTrue(max(dimSums(out[, , landShares], 3)) <= 1 + 10^-5,
                 "land shares sum to <= 1 + 10^-5") # if < 1 the rest is likely water
  toolExpectTrue(all(out[, -1, c("primf", "primn")] <= setYears(out[, -nyears(out), c("primf", "primn")],
                                                                getYears(out[, -1, ]))),
                 "primf and primn are never expanding", falseStatus = "warn")

  return(list(x = out,
              isocountries = FALSE,
              unit = "1",
              min = 0,
              max = 1,
              description = paste("land use data downscaled to 0.25 degree")))
}
