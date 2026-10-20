#' fullSCENARIOMIP
#'
#' Run the pipeline to generate harmonized and downscaled data to report for ScenarioMIP.
#' LUH3 is used as historical reference dataset for harmonization and downscaling.
#' Write .nc files, print full report on consistency checks and write it to report.log.
#'
#' @param rev revision number of the data. If not provided the current date will be used instead.
#' When called via madrat::retrieveData rev will be converted to numeric_version.
#' @inheritParams calcLandInput
#' @param scenario scenario name to be included in filenames
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @param yearsSubset remove years from the returned data which are not in yearsSubset
#' @param harmonization name of harmonization method, see \code{\link{toolGetHarmonizer}}
#' @param downscaling name of downscaling method: "proportional" (\code{luscale::interpolate2}) or "transitions"
#' (\code{\link{toolDownscaleTransitions}})
#' @param compression compression level of the resulting .nc files, possible values are integers from 1-9,
#' 1 = fastest, 9 = best compression
#' @param progress boolean defining whether progress should be printed
#'
#' @author Pascal Sauer
fullSCENARIOMIP <- function(rev = numeric_version("0"), input = "magpie", scenario = "",
                            harmonizationPeriod = c(2025, 2050),
                            yearsSubset = 1995:2100,
                            harmonization = "fadeForest", downscaling = "proportional",
                            compression = 2, progress = TRUE) {
  revision <- if (identical(rev, numeric_version("0"))) format(Sys.time(), "%Y-%m-%d") else rev

  fileSuffix <- paste0("_input4MIPs_landState_ScenarioMIP_",
                       scenario, if (scenario != "") "-",
                       revision, "_gn_", min(yearsSubset), "-", max(yearsSubset), ".nc")

  writeArgs <- list(compression = compression, missval = 1e20, progress = progress,
                    gridDefinition = c(-179.875, 179.875, -89.875, 89.875, 0.25))

  metadataArgs <- list(revision = revision, missingValue = 1e20, resolution = 0.25,
                       compression = compression, harmonizationPeriod = harmonizationPeriod,
                       activityId = "ScenarioMIP",
                       references = paste0("https://github.com/pik-piam/mrdownscale and ",
                                           "https://wcrp-cmip.org/mips/scenariomip/"),
                       targetMIP = "ScenarioMIP",
                       ncTitle = paste0("REMIND-MAgPIE Land-Use Data Harmonized ",
                                        "and Downscaled using LUH3 historic as reference"),
                       referenceDataset = paste0("LUH3 historic from https://aims2.llnl.gov/search/input4mips/ ",
                                                 "(institution_id = 'UofMD' and mip_era = 'CMIP7')"),
                       furtherInfoUrl = "NA")

  ncFile <- paste0("multiple-states", fileSuffix)
  calcOutput("StatesNC", outputFormat = "ScenarioMIP", input = input,
             harmonizationPeriod = harmonizationPeriod,
             yearsSubset = yearsSubset,
             harmonization = harmonization, downscaling = downscaling,
             aggregate = FALSE, file = ncFile, writeArgs = writeArgs)
  do.call(toolAddMetadataNC, c(ncFile = ncFile, metadataArgs))

  ncFile <- paste0("multiple-management", fileSuffix)
  calcOutput("ManagementNC", outputFormat = "ScenarioMIP", input = input,
             harmonizationPeriod = harmonizationPeriod,
             yearsSubset = yearsSubset,
             harmonization = harmonization, downscaling = downscaling,
             aggregate = FALSE, file = ncFile, writeArgs = writeArgs)
  do.call(toolAddMetadataNC, c(ncFile = ncFile, metadataArgs))

  if (input == "magpie") {
    ncFile <- paste0("multiple-transitions", fileSuffix)
    calcOutput("TransitionsNC", outputFormat = "ScenarioMIP", input = input,
               harmonizationPeriod = harmonizationPeriod,
               yearsSubset = yearsSubset,
               harmonization = harmonization, downscaling = downscaling,
               aggregate = FALSE, file = ncFile, writeArgs = writeArgs)
    do.call(toolAddMetadataNC, c(ncFile = ncFile, metadataArgs))
  }

  toolWriteMadratLog(logPath = "consistencyCheck.log")
}
