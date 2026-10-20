#' calcAgeClassesReport
#'
#' Age classes on the grid (\code{\link{calcAgeClassesHighRes}}) in the format of an output.
#'
#' @param outputFormat "downscaledmagpie": for each group the share of its area in each age class (NA
#' where the grid cell has no area of the group), and the diagnostics
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
    out[[group]] <- ages[[group]] / dimSums(ages[[group]], dim = 3)
    out[[group]][is.nan(out[[group]])] <- NA
  }
  out[["diagnostics"]] <- ages[["diagnostics"]]
  return(list(x = out, class = "list", unit = "1", cache = FALSE,
              description = "Share of the area of each group in each age class by grid cell"))
}
