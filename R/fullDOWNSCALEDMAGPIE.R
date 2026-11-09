#' fullDOWNSCALEDMAGPIE
#'
#' Run the pipeline to generate harmonized and downscaled MAgPIE data using
#' landuseinit as a reference dataset. Write output in the format of
#' avl_land_t_0.5.mz, full report on consistency checks is printed and
#' written to report.log.
#'
#' @param rev revision number of the data. If not provided the current date will be used instead.
#' When called via madrat::retrieveData rev will be converted to numeric_version.
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @param target Name of dataset to be used as harmonization target and downscaling reference
#' @param harmonization name of harmonization method, see \code{\link{toolGetHarmonizer}}
#' @param downscaling name of downscaling method: "proportional" (\code{luscale::interpolate2}) or "transitions"
#' (\code{\link{toolDownscaleTransitions}})
#'
#' @author Pascal Sauer
fullDOWNSCALEDMAGPIE <- function(rev = numeric_version("0"),
                                 harmonizationPeriod = c(2015, 2050), target = "landuseinit",
                                 downscaling = "proportional", harmonization = "fade") {
  if (target != "landuseinit") {
    stop("currently only target = landuseinit is supported")
  }

  # LandReport with outputFormat = downscaledmagpie assumes target = landuseinit
  calcOutput("LandReport", outputFormat = "downscaledmagpie", input = "magpie",
             harmonizationPeriod = harmonizationPeriod, yearsSubset = seq(1995, 2100, 5),
             harmonization = harmonization, downscaling = downscaling,
             aggregate = FALSE, file = "cell.land_0.5.mz")

  toolWriteMadratLog(logPath = "consistencyCheck.log")
}
