#' calcAgeClassesReport
#'
#' Age classes on the grid (\code{\link{calcAgeClassesHighRes}}) in the format of an output.
#'
#' @param outputFormat "downscaledmagpie": for each group the share of its area in each age class (NA
#' where the grid cell has less than 1e-9 Mha of the group, where the composition is undefined), and the
#' diagnostics
#' @inheritParams calcAgeClassesHighRes
#' @return List with a magclass object of shares per group and \code{diagnostics}
#' @author Patrick v. Jeetze
#' @seealso \code{\link{fullDOWNSCALEDMAGPIE}}
calcAgeClassesReport <- function(outputFormat, input, harmonization = "none") {
  if (outputFormat != "downscaledmagpie") {
    stop("Can only report for outputFormat = downscaledmagpie")
  }
  ages <- calcOutput("AgeClassesHighRes", input = input, harmonization = harmonization, aggregate = FALSE)

  out <- list()
  for (group in setdiff(names(ages), "diagnostics")) {
    area <- dimSums(ages[[group]], dim = 3)
    out[[group]] <- ages[[group]] / area
    # below 1e-9 Mha of the group (its area repeated for each age class), the composition is rounding noise
    out[[group]][ages[[group]] * 0 + area < 1e-9] <- NA
  }
  out[["diagnostics"]] <- ages[["diagnostics"]]
  return(list(x = out, class = "list", unit = "1", cache = FALSE,
              description = "Share of the area of each group in each age class by grid cell"))
}
