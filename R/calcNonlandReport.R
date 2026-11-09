#' calcNonlandReport
#'
#' Convert the downscaled nonland data to the format required by the given project.
#'
#' @param outputFormat options: ESM, ScenarioMIP
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @param yearsSubset vector of years to keep in the output dataset
#' @param harmonization name of harmonization method, see \code{\link{toolGetHarmonizer}}
#' @param downscaling name of downscaling method: "proportional" (\code{luscale::interpolate2}) or "transitions"
#' (\code{\link{toolDownscaleTransitions}})
#' @return nonland data
#'
#' @examples
#' \dontrun{
#'   calcOutput("NonlandReport", outputFormat = "ESM",
#'              harmonizationPeriod = c(2015, 2050), yearsSubset = 2015:2100,
#'              harmonization = "fade", downscaling = "proportional")
#' }
#' @author Pascal Sauer
calcNonlandReport <- function(outputFormat, harmonizationPeriod, yearsSubset, harmonization, downscaling) {
  if (outputFormat %in% c("ESM", "ScenarioMIP")) {
    input <- "magpie"
    if (outputFormat == "ESM") {
      target <- "luh2mod"
      cellAreaKm2 <- readSource("LUH2v2h", subtype = "cellArea", convert = FALSE)
    } else if (outputFormat == "ScenarioMIP") {
      target <- "luh3"
      cellAreaKm2 <- readSource("LUH3", subtype = "cellArea", convert = FALSE)
    }
    x <- calcOutput("NonlandHighRes", input = input, target = target,
                    harmonizationPeriod = harmonizationPeriod, yearsSubset = yearsSubset,
                    harmonization = harmonization, downscaling = downscaling, aggregate = FALSE)

    if (outputFormat == "ScenarioMIP") {
      # combine secyf + secmf into secdf
      secyfWha <- "wood_harvest_area.secyf"
      secmfWha <- "wood_harvest_area.secmf"
      x[, , secyfWha] <- dimSums(x[, , c(secyfWha, secmfWha)], dim = 3)
      x <- x[, , secmfWha, invert = TRUE]
      getNames(x) <- sub(secyfWha, "wood_harvest_area.secdf", getNames(x))

      x[, , "bioh.secyf"] <- dimSums(x[, , c("bioh.secyf", "bioh.secmf")], dim = 3)
      x <- x[, , "bioh.secmf", invert = TRUE]
      getNames(x) <- sub("bioh.secyf", "bioh.secdf", getNames(x))
    }

    cellAreaKm2 <- as.magpie(cellAreaKm2)
    stopifnot(getItems(x, 1) %in% getItems(cellAreaKm2, 1))
    cellAreaKm2 <- collapseDim(cellAreaKm2[getItems(x, 1), , ], 3)
    # convert from km2 to Mha
    cellAreaMha <- cellAreaKm2 / 10000

    # convert from Mha to shares
    harv <- x[, , "wood_harvest_area"] / cellAreaMha
    getNames(harv) <- sub("wood_harvest_area", "harv", getNames(harv))

    woodTypeShares <- x[, , "harvest_weight_type"]

    # get rndwd/fulwd shares per country to replace NAs
    coords <- getCoords(woodTypeShares)
    mapping <- calcOutput("ResolutionMapping", input = input, target = target, aggregate = FALSE)
    mapping <- mapping[, c("x", "y", "country")]
    merged <- merge(coords, mapping, sort = FALSE)
    stopifnot(identical(merged[, c("x", "y")], coords))
    getItems(woodTypeShares, 1, raw = TRUE) <- paste0(getItems(woodTypeShares, 1), ".", merged$country)
    names(dimnames(woodTypeShares))[1] <- "x.y.country"
    countryTotals <- toolAggregate(woodTypeShares, to = "country")
    countryWoodTypeShares <- countryTotals / dimSums(countryTotals, 3)
    countryWoodTypeShares[is.na(countryWoodTypeShares)] <- 0.5 # replace NAs with share 0.5, so sum is still 1

    total <- dimSums(woodTypeShares, 3) # 1269815 cells with total == 0
    woodTypeShares <- woodTypeShares / total # NAs introduced by cells with total == 0
    stopifnot(sum(is.na(woodTypeShares)) == 2 * sum(total == 0))
    fillValuesCountry <- countryWoodTypeShares[getItems(woodTypeShares, "country", full = TRUE), , ]
    getItems(fillValuesCountry, 1, raw = TRUE) <- getItems(woodTypeShares, 1)
    names(dimnames(fillValuesCountry))[[1]] <- "x.y.country"
    woodTypeShares[is.na(woodTypeShares)] <- fillValuesCountry[is.na(woodTypeShares)]
    getNames(woodTypeShares) <- c("harvest_weight_type.rndwd", "harvest_weight_type.fulwd")
    woodTypeShares[, , "rndwd"] <- 1 - woodTypeShares[, , "fulwd"]
    woodTypeShares <- collapseDim(woodTypeShares, 1.3)

    fertl <- x[, , "fertilizer"]
    getNames(fertl) <- sub("fertilizer", "fertl", getNames(fertl))
    out <- mbind(fertl, harv, woodTypeShares, x[, , "bioh"])
    getNames(out) <- sub("\\.", "_", getNames(out))
    getNames(out) <- sub("^bioh_(.+)", "\\1_bioh", getNames(out))
    getNames(out) <- sub("^harv_(.+)", "\\1_harv", getNames(out))
    getNames(out) <- sub("harvest_weight_type_", "", getNames(out))
    names(dimnames(out))[3] <- "data"

    shares <- grep("(harv|rndwd|fulwd)$", getNames(out), value = TRUE)
    if (max(out[, , shares]) > 1.0001) {
      toolStatusMessage("warn", paste0("Some shares are > 1 (max: ", max(out[, , shares]), "), setting those to 1"))
    }
    out[, , shares][out[, , shares] > 1] <- 1

    toolExpectTrue(min(out) >= 0, "All values are >= 0")
    toolExpectTrue(max(out[, , shares]) <= 1, "All shares are <= 1")
    toolExpectLessDiff(dimSums(out[, , c("rndwd", "fulwd")], 3), 1, 10^-5, "shares sum to 1 (rndwd + fulwd)")
    toolExpectTrue(!any(is.na(out)), "No NAs in output")

    return(list(x = out,
                isocountries = FALSE,
                min = 0,
                unit = "rndwd & fulwd: 1; bioh: kg C yr-1; harv: 1; fertl: kg ha-1 yr-1",
                description = "Downscaled nonland data report for use in ESMs"))
  } else {
    stop("Can only report for outputFormat = ESM/ScenarioMIP")
  }
}
