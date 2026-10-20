#' fullESM
#'
#' Run the pipeline to generate harmonized and downscaled data to report for RESCUE
#' and other projects where ESM compatible land use inputs are required.
#' Write .nc files, print full report on consistency checks and write it to report.log.
#'
#' @param rev revision number of the data. If not provided the current date will be used instead.
#' When called via madrat::retrieveData rev will be converted to numeric_version.
#' @inheritParams calcLandInput
#' @param fileNamePart scenario nc files will be called
#' multiple-[type]_[fileNamePart]-[revision]_gn_[min(yearsSubset)]-[max(yearsSubset)].nc
#' where type is one of states/management/transitions
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
#' @author Pascal Sauer, Jan Philipp Dietrich
fullESM <- function(rev = numeric_version("0"), input = "magpie", fileNamePart = "",
                    harmonizationPeriod = c(2015, 2050), yearsSubset = 2015:2100,
                    harmonization = "fade", downscaling = "proportional",
                    compression = 2, progress = TRUE) {
  revision <- if (identical(rev, numeric_version("0"))) format(Sys.time(), "%Y-%m-%d") else rev

  fileSuffix <- paste0(fileNamePart, if (fileNamePart != "") "-",
                       revision, "_gn_", min(yearsSubset), "-", max(yearsSubset), ".nc")

  writeArgs <- list(compression = compression, missval = 1e20, progress = progress,
                    gridDefinition = c(-179.875, 179.875, -89.875, 89.875, 0.25))

  metadataArgs <- list(revision = revision, missingValue = 1e20, resolution = 0.25,
                       compression = compression, harmonizationPeriod = harmonizationPeriod,
                       activityId = "RESCUE",
                       references = paste0("https://github.com/pik-piam/mrdownscale and https://rescue-climate.eu"),
                       targetMIP = "RESCUE",
                       ncTitle = "MAgPIE Land-Use Data Harmonized and Downscaled using LUH2 v2h as reference",
                       referenceDataset = "LUH2 v2h Release (10/14/16) from https://luh.umd.edu/data.shtml",
                       furtherInfoUrl = paste0("https://github.com/pik-piam/mrdownscale/blob/",
                                               "f358ccc1902da769e49c5d52e1085db6b8c797b3/changelog.md"))

  ncFile <- paste0("multiple-states_", fileSuffix)
  calcOutput("StatesNC", outputFormat = "ESM", input = input,
             harmonizationPeriod = harmonizationPeriod,
             yearsSubset = yearsSubset,
             harmonization = harmonization, downscaling = downscaling,
             aggregate = FALSE, file = ncFile, writeArgs = writeArgs)
  do.call(toolAddMetadataNC, c(ncFile = ncFile, metadataArgs))

  ncFile <- paste0("multiple-management_", fileSuffix)
  calcOutput("ManagementNC", outputFormat = "ESM", input = input,
             harmonizationPeriod = harmonizationPeriod,
             yearsSubset = yearsSubset,
             harmonization = harmonization, downscaling = downscaling,
             aggregate = FALSE, file = ncFile, writeArgs = writeArgs)
  do.call(toolAddMetadataNC, c(ncFile = ncFile, metadataArgs))

  if (input == "magpie") {
    ncFile <- paste0("multiple-transitions_", fileSuffix)
    calcOutput("TransitionsNC", outputFormat = "ESM", input = input,
               harmonizationPeriod = harmonizationPeriod,
               yearsSubset = yearsSubset,
               harmonization = harmonization, downscaling = downscaling,
               aggregate = FALSE, file = ncFile, writeArgs = writeArgs)
    do.call(toolAddMetadataNC, c(ncFile = ncFile, metadataArgs))
  }

  toolWriteMadratLog(logPath = "consistencyCheck.log")
}
