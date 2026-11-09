#' calcStatesNC
#'
#' Prepare data to be written as LUH-style states.nc file. Call this via calcOutput
#' in a full function, and set calcOutput's file argument to a .nc file path.
#'
#' @param outputFormat options: ESM, ScenarioMIP
#' @inheritParams calcLandInput
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @param yearsSubset remove years from the returned data which are not in yearsSubset
#' @param harmonization name of harmonization method, see \code{\link{toolGetHarmonizer}}
#' @param downscaling name of downscaling method: "proportional" (\code{luscale::interpolate2}) or "transitions"
#' (\code{\link{toolDownscaleTransitions}})
#' @return data prepared to be written as a LUH-style states.nc file
#' @author Pascal Sauer, Jan Philipp Dietrich
calcStatesNC <- function(outputFormat, input, harmonizationPeriod, yearsSubset, harmonization, downscaling) {
  x <- calcOutput("LandReport", outputFormat = outputFormat, input = input,
                  harmonizationPeriod = harmonizationPeriod, yearsSubset = yearsSubset,
                  harmonization = harmonization, downscaling = downscaling, aggregate = FALSE)

  statesVariables <- c("c3ann", "c3nfx", "c3per", "c4ann", "c4per", "pastr",
                       "primf", "primn", "range", "secdf", "secdn", "urban")
  if (outputFormat == "ScenarioMIP") {
    statesVariables <- c(statesVariables, "pltns")
  }
  x <- x[, , statesVariables]

  # account for unit "years since 1970-01-01 0:0:0"
  x <- setYears(x, getYears(x, as.integer = TRUE) - 1970)

  if (outputFormat == "ScenarioMIP") {
    expectedVariables <- c("primf", "primn", "secdf", "secdn", "pastr", "range", "urban",
                           "c3ann", "c3per", "c4ann", "c4per", "c3nfx", "pltns")
    toolExpectTrue(setequal(getItems(x, 3), expectedVariables), "variable names are as expected")
  }

  return(list(x = x,
              isocountries = FALSE,
              unit = "1",
              min = 0,
              max = 1,
              cache = FALSE,
              description = "data prepared to be written as a LUH-style states.nc file"))
}
