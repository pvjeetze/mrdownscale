#' calcNonlandHighRes
#'
#' Calculate a high resolution dataset from the low resolution input dataset
#' and high resolution data.
#'
#' Wood harvest area is disaggregated using the maximum possible harvest per
#' year which is based on the downscaled land data. Bioh
#' is disaggregated using the just disaggregated wood harvest area as weight.
#' Fertilizer is in kg ha-1 yr-1, so we simply use the low res/region value for each
#' cell corresponding to that region. Harvest weight type is disaggregated
#' using the nonland target data in the first year of
#' the harmonization period as weight.
#'
#' @inheritParams calcNonlandInput
#' @param target name of a target dataset, see \code{\link{calcLandTarget}}
#' for available target datasets
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @param yearsSubset vector of years to keep in the output dataset
#' @param harmonization name of harmonization method, see \code{\link{toolGetHarmonizer}}
#' @param downscaling name of downscaling method: "proportional" (\code{luscale::interpolate2}) or "transitions"
#' (\code{\link{toolDownscaleTransitions}})
#' @return downscaled nonland data
#'
#' @examples
#' \dontrun{
#'   calcOutput("NonlandHighRes", input = "magpie", target = "luh2mod",
#'              harmonizationPeriod = c(2015, 2050), yearsSubset = 2015:2100,
#'              harmonization = "fade", downscaling = "proportional")
#' }
#' @author Pascal Sauer
calcNonlandHighRes <- function(input, target, harmonizationPeriod, yearsSubset, harmonization, downscaling) {
  hp <- harmonizationPeriod
  x <- calcOutput("NonlandHarmonized", input = input, target = target,
                  harmonizationPeriod = hp, harmonization = harmonization, aggregate = FALSE)
  x <- x[, getYears(x, as.integer = TRUE) %in% yearsSubset, ]

  nonlandTarget <- calcOutput("NonlandTarget", target = target, endOfHistory = hp[1], aggregate = FALSE)
  nonlandTarget <- as.magpie(nonlandTarget[[terra::time(nonlandTarget) %in% yearsSubset]])
  getItems(nonlandTarget, 3, raw = TRUE) <- sub("^(.+?)_(.+)$", "\\2.\\1", getItems(nonlandTarget, 3))
  names(dimnames(nonlandTarget))[3] <- "category.data"

  futureYears <- setdiff(getYears(x, TRUE), getYears(nonlandTarget, TRUE))

  resmap <- calcOutput("ResolutionMapping", input = input, target = target, aggregate = FALSE)

  landHighRes <- calcOutput("LandHighRes", input = input, target = target,
                            harmonizationPeriod = hp, yearsSubset = yearsSubset,
                            harmonization = harmonization, downscaling = downscaling, aggregate = FALSE)

  message("downscaling wood harvest area...")
  harvestArea <- x[, futureYears, "wood_harvest_area"]
  harvestAreaWeight <- toolMaxHarvestPerYear(landHighRes)[, futureYears, ]
  harvestAreaDownscaled <- toolAggregate(harvestArea, resmap, weight = harvestAreaWeight,
                                         from = "lowRes", to = "cell", dim = 1, zeroWeight = "fix")

  message("downscaling wood harvest mass...")
  bioh <- x[, futureYears, "bioh"]
  weightBioh <- harvestAreaDownscaled
  getItems(weightBioh, 3.1) <- sub("wood_harvest_area$", "bioh", getItems(weightBioh, 3.1))
  stopifnot(getItems(weightBioh, 3) == getItems(bioh, 3))
  biohDownscaled <- toolAggregate(bioh, resmap, weight = weightBioh, from = "lowRes", to = "cell",
                                  dim = 1, zeroWeight = "fix")

  # fertilizer is in kg ha-1 yr-1 already, use low res/region value for each cell corresponding to that region
  fertilizerDownscaled <- setCells(x[resmap$lowRes, futureYears, "fertilizer"], resmap$cell)

  message("downscaling wood harvest weight type...")
  harvestType <- x[, futureYears, "harvest_weight_type"]
  weightHarvestType <- collapseDim(nonlandTarget[, nyears(nonlandTarget), "harvest_weight_type"])
  harvestTypeDownscaled <- toolAggregate(harvestType, resmap, weight = weightHarvestType,
                                         from = "lowRes", to = "cell", dim = 1, zeroWeight = "fix")


  out <- mbind(nonlandTarget,
               mbind(harvestAreaDownscaled, biohDownscaled, fertilizerDownscaled, harvestTypeDownscaled))

  # checks
  fertilizerInput <- calcOutput("NonlandInputRecategorized", input = input, target = target, aggregate = FALSE)

  inSum <- dimSums(x[, , "fertilizer", invert = TRUE], dim = 1)
  outSum <- dimSums(out[, , "fertilizer", invert = TRUE], dim = 1)
  stopifnot(identical(getYears(inSum), getYears(outSum)),
            setequal(getItems(inSum, 3), getItems(outSum, 3)))
  toolExpectLessDiff(inSum, outSum, 10^-3,
                     "No significant global sum difference per category before and after downscaling")

  # for years after harmonization make sure that total global fertilizer matches input
  years <- getYears(out, TRUE)
  years <- years[years >= hp[2]]
  fertilizerInput <- dimSums(fertilizerInput[, years, "fertilizer"], 1)
  fertilizerOutput <- toolFertilizerTg(out[, years, "fertilizer"], landHighRes[, years, ])
  toolExpectLessDiff(fertilizerInput, dimSums(fertilizerOutput, 1), 10^-5,
                     "Total global fertilizer after harmonization period matches input data")
  toolCheckFertilizer(out[, , "fertilizer"])

  toolExpectTrue(setequal(getItems(out, dim = 3), getItems(x, dim = 3)),
                 "Nonland categories remain unchanged")
  toolExpectTrue(min(out) >= 0, "All values are >= 0")

  toolCheckWoodHarvestArea(out[, , "wood_harvest_area"], landHighRes, hp[1])

  return(list(x = out,
              min = 0,
              isocountries = FALSE,
              unit = "harvest_weight & bioh: kg C yr-1; harvest_area: Mha yr-1; fertilizer: kg ha-1 yr-1",
              description = "Downscaled nonland data"))
}
