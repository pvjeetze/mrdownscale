#' calcManagementNC
#'
#' Prepare data to be written as a LUH-style management.nc file.  Call this via calcOutput
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
#' @return data prepared to be written as a LUH-style management.nc file
#' @author Pascal Sauer, Jan Philipp Dietrich
calcManagementNC <- function(outputFormat, input, harmonizationPeriod, yearsSubset, harmonization, downscaling) {
  x <- calcOutput("LandReport", outputFormat = outputFormat, input = input,
                  harmonizationPeriod = harmonizationPeriod, yearsSubset = yearsSubset,
                  harmonization = harmonization, downscaling = downscaling, aggregate = FALSE)

  if (outputFormat == "ESM") {
    landManagementVariables <- c("crpbf_c3ann", "crpbf_c3nfx", "crpbf_c3per", "crpbf_c4ann", "crpbf_c4per",
                                 "crpbf2_c3per", "crpbf2_c4per",
                                 "irrig_c3ann", "irrig_c3nfx", "irrig_c3per", "irrig_c4ann", "irrig_c4per",
                                 "manaf")
  } else if (outputFormat == "ScenarioMIP") {
    landManagementVariables <- c("cpbf1_c3ann", "cpbf1_c3nfx", "cpbf1_c3per", "cpbf1_c4ann", "cpbf1_c4per",
                                 "cpbf2_c3per", "cpbf2_c4per",
                                 "irrig_c3ann", "irrig_c3nfx", "irrig_c3per", "irrig_c4ann", "irrig_c4per")
  }
  x <- x[, , intersect(landManagementVariables, getItems(x, 3))]

  if (input == "magpie") {
    nonland <- calcOutput("NonlandReport", outputFormat = outputFormat,
                          harmonizationPeriod = harmonizationPeriod,
                          yearsSubset = yearsSubset,
                          harmonization = harmonization, downscaling = downscaling,
                          aggregate = FALSE)
    nonlandManagementVariables <- c("fertl_c3nfx", "fertl_c3per", "fertl_c3ann", "fertl_c4ann", "fertl_c4per",
                                    "rndwd", "fulwd")
    nonland <- nonland[, , nonlandManagementVariables]

    if (outputFormat == "ScenarioMIP") {
      nonland <- add_columns(nonland, "pltns_wdprd", fill = 1)
      nonland <- add_columns(nonland, "pltns_bfuel", fill = 0)
    }
    x <- mbind(x, nonland)
  }

  # account for the time unit written into the nc file by terra: "years since 1970-01-01 0:0:0"
  x <- setYears(x, getYears(x, as.integer = TRUE) - 1970)

  unit <- "1"
  if (outputFormat == "ScenarioMIP") {
    expectedVariables <- c("irrig_c3ann", "irrig_c3per", "irrig_c4ann", "irrig_c4per", "irrig_c3nfx",
                           "cpbf1_c3ann", "cpbf1_c4ann", "cpbf1_c3per", "cpbf1_c4per", "cpbf1_c3nfx",
                           "cpbf2_c3per", "cpbf2_c4per")
    if (input == "magpie") {
      expectedVariables <- c(expectedVariables,
                             "fertl_c3ann", "fertl_c4ann", "fertl_c3per", "fertl_c4per", "fertl_c3nfx",
                             "rndwd", "fulwd",
                             "pltns_wdprd", "pltns_bfuel")
      unit <- "1, except fertl: kg ha-1 yr-1"
    }

    toolExpectTrue(setequal(getItems(x, 3), expectedVariables), "variable names are as expected")
  }

  return(list(x = x,
              isocountries = FALSE,
              unit = unit,
              min = 0,
              cache = FALSE,
              description = "data prepared to be written as a LUH-style management.nc file"))
}
