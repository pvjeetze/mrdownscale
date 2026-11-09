#' calcTransitionsNC
#'
#' Prepared data to be written as a LUH-style transitions.nc file
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
#' @return data prepared to be written as a LUH-style transitions.nc file
#' @author Pascal Sauer, Jan Philipp Dietrich
calcTransitionsNC <- function(outputFormat, input, harmonizationPeriod, yearsSubset, harmonization, downscaling) {
  nonland <- calcOutput("NonlandReport", outputFormat = outputFormat,
                        harmonizationPeriod = harmonizationPeriod, yearsSubset = yearsSubset,
                        harmonization = harmonization, downscaling = downscaling, aggregate = FALSE)

  nonland <- nonland[, , grep("_(bioh|harv)$", getItems(nonland, 3), value = TRUE)]

  x <- nonland[, getYears(nonland, as.integer = TRUE) %in% yearsSubset, ]

  if (outputFormat == "ESM") {
    transitions <- calcOutput("LandTransitions", outputFormat = outputFormat,
                              harmonizationPeriod = harmonizationPeriod, yearsSubset = yearsSubset,
                              harmonization = harmonization, downscaling = downscaling,
                              input = input, gross = TRUE, aggregate = FALSE)
    getItems(transitions, raw = TRUE, dim = 3) <- sub("\\.", "_to_", getItems(transitions, dim = 3))
    getSets(transitions, fulldim = FALSE)[3] <- "transitions"
    # we use to-semantics for transitions (value for 1994 describes what happens from 1993 to 1994)
    # by subtracting 1 we get from-semantics (value for 1994 describes what happens from 1994 to 1995)
    # which is what LUH uses
    getYears(transitions) <- getYears(transitions, as.integer = TRUE) - 1
    transitions <- transitions[, getYears(transitions, as.integer = TRUE) %in% yearsSubset, ]
    x <- mbind(x, transitions)
  }

  # account for unit "years since 1970-01-01 0:0:0"
  x <- setYears(x, getYears(x, as.integer = TRUE) - 1970)

  if (outputFormat == "ScenarioMIP") {
    expectedVariables <- c("primf_harv", "secdf_harv", "primn_harv", "secnf_harv", "pltns_harv",
                           "primf_bioh", "secdf_bioh", "primn_bioh", "secnf_bioh", "pltns_bioh")
    toolExpectTrue(setequal(getItems(x, 3), expectedVariables), "variable names are as expected")
  }

  return(list(x = x,
              isocountries = FALSE,
              unit = "1 or kg C yr-1",
              min = 0,
              cache = FALSE,
              description = "data prepared to be written as a LUH-style transitions.nc file"))
}
