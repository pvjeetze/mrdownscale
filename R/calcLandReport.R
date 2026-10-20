#' calcLandReport
#'
#' Convert the downscaled land use data to the format required by the given project.
#'
#' @param outputFormat format in which the outputs should be prepared. Options: ESM, ScenarioMIP,
#' downscaledmagpie (land by land type, as avl_land_t_0.5.mz), and with \code{downscaling = "transitions"} and
#' \code{harmonization = "none"} downscaledmagpieSplit (land by land type, sub-type and water supply, names
#' "land.type.water": crop area of food and feed crops kfo and of bioenergy crops kbe, rainfed and
#' irrigated, fallow land and tree cover, and the sub-types of the other land types with water "total")
#' and downscaledmagpieCroparea (crop area by crop and water supply as share of the grid cell area, in
#' the years of the model)
#' @inheritParams calcLandInput
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @param yearsSubset vector of years to keep in the output dataset
#' @inheritParams calcLandHighRes
#' @return land use data
#' @author Pascal Sauer, Patrick v. Jeetze
calcLandReport <- function(outputFormat, input, harmonizationPeriod, yearsSubset, harmonization, downscaling) {
  if (outputFormat == "ESM") {
    native <- calcOutput("LandHighRes", input = input, target = "luh2mod",
                         harmonizationPeriod = harmonizationPeriod, yearsSubset = yearsSubset,
                         downscaling = downscaling, harmonization = harmonization,
                         aggregate = FALSE)
    cellArea <- readSource("LUH2v2h", subtype = "cellArea", convert = FALSE)
    cellArea <- collapseDim(as.magpie(cellArea), 3)

    cropData <- toolCropData(native, cellArea)
    getNames(cropData) <- sub("^cpbf1_", "crpbf_", getNames(cropData))
    getNames(cropData) <- sub("^cpbf2_", "crpbf2_", getNames(cropData))

    totalSecondaryForest <- dimSums(native[, , c("pltns", "secdf")], dim = 3)
    # calculate manaf (managed forest) = pltns share of total secondary forest
    manaf <- native[, , "pltns"] / totalSecondaryForest
    manaf[totalSecondaryForest == 0] <- 0 # replace NAs introduced by 0 / 0
    getNames(manaf) <- "manaf"

    nonCropData <- native[, , c("primf", "primn", "secdf", "secdn", "urban", "pastr", "range")]
    nonCropData[, , "secdf"] <- totalSecondaryForest
    nonCropData <- nonCropData * 10000 / cellArea[getItems(cellArea, 1) %in% getItems(nonCropData, 1), , ]

    out <- mbind(cropData, nonCropData, manaf)

    if (max(out) > 1.0001) {
      toolStatusMessage("warn", paste0("Some shares are > 1 (max: ", max(out), "), setting those to 1"))
    }
    out[out > 1] <- 1

    toolExpectTrue(min(out) >= 0, "All values are >= 0")
    toolExpectTrue(max(out) <= 1, "All shares are <= 1")
    landShares <- c("c3ann", "c3nfx",  "c3per", "c4ann", "c4per", "primf",
                    "primn", "secdf", "secdn", "urban", "pastr", "range")
    toolExpectTrue(max(dimSums(out[, , landShares], 3)) <= 1 + 10^-5,
                   "land shares sum to <= 1 + 10^-5") # if < 1 the rest is likely water
    toolExpectTrue(all(out[, -1, c("primf", "primn")] <= setYears(out[, -nyears(out), c("primf", "primn")],
                                                                  getYears(out[, -1, ]))),
                   "primf and primn are never expanding", falseStatus = "warn")

    return(list(x = out,
                isocountries = FALSE,
                unit = "1",
                min = 0,
                max = 1,
                description = paste("land use data downscaled to LUH2 resolution")))
  } else if (outputFormat == "ScenarioMIP") {
    report <- calcOutput("LandReportScenarioMIP",
                         input = input,
                         harmonizationPeriod = harmonizationPeriod,
                         yearsSubset = yearsSubset,
                         harmonization = harmonization, downscaling = downscaling,
                         supplementary = TRUE, aggregate = FALSE)
    return(list(x = report$x,
                isocountries = report$isocountries,
                unit = report$unit,
                min = report$min,
                max = report$max,
                description = report$description))
  } else if (outputFormat == "downscaledmagpie") {
    target <- "landuseinit"
    # along the transitions of the run: all years of the model, the start as y1985 and the grid cells of
    # the run included (see calcLandHighRes)
    transitions <- downscaling == "transitions" && harmonization == "none"
    yearsSubset <- seq(1995, 2100, 5)
    if (transitions) {
      yearsSubset <- NULL
    }
    x <- calcOutput("LandHighRes", input = input, target = target,
                    harmonizationPeriod = harmonizationPeriod, yearsSubset = yearsSubset,
                    downscaling = downscaling, harmonization = harmonization, aggregate = FALSE)

    if (!transitions) {
      # add iso
      mapping <- calcOutput("ResolutionMapping", input = input, target = target, aggregate = FALSE)
      stopifnot(all(getItems(x, 1) %in% mapping$cell))
      mapping <- mapping[match(getItems(x, 1), mapping$cell), c("cell", "country")]
      getItems(x, 1, raw = TRUE) <- paste0(mapping$cell, ".", mapping$country)

      # rename primf/secdf to primforest/secdforest
      getItems(x, 3) <- sub("primf", "primforest", getItems(x, 3))
      getItems(x, 3) <- sub("secdf", "secdforest", getItems(x, 3))

      # original unchanged LanduseInit 1995 -> 1985 in output
      x <- mbind(setYears(readSource("LanduseInit")[, 1995, ],
                          1985),
                 x)
    }

    names(dimnames(x)) <- c("x.y.iso", "year", "data")

    return(list(x = x,
                isocountries = FALSE,
                unit = "Mha",
                min = 0,
                max = 1,
                description = paste0("MAgPIE data harmonized and downscaled using landuseinit as reference ",
                                     "for further processing in magpie4")))
  } else if (outputFormat %in% c("downscaledmagpieSplit", "downscaledmagpieCroparea")) {
    if (!(downscaling == "transitions" && harmonization == "none")) {
      stop("outputFormat = \"", outputFormat, "\" needs downscaling = \"transitions\" and harmonization = \"none\"")
    }
    land <- calcOutput("LandTransitionsHighRes", input = input, harmonization = harmonization, aggregate = FALSE)
    landHr <- land[["landHr"]]
    cropSplitHr <- land[["splitHr"]][["crop"]]
    cropAreaHr <- land[["splitHr"]][["crop_area"]]
    rm(land)

    if (outputFormat == "downscaledmagpieCroparea") {
      # crop area by crop and water supply as share of the grid cell area, in the years of the model
      x <- cropAreaHr / setNames(dimSums(landHr, dim = 3), NULL)
      x <- x[, getYears(x)[-1], ]
      getSets(x, fulldim = FALSE) <- c("x.y.iso", "year", "kcr.w")
      return(list(x = x, isocountries = FALSE, unit = "1", min = 0,
                  description = "MAgPIE crop area by crop and water supply as share of the grid cell area"))
    }

    # crop area of food and feed crops (kfo) and of bioenergy crops (kbe, betr and begr), each rainfed and
    # irrigated, fallow land and tree cover; the other land types by sub-type
    kbe <- c("betr", "begr")
    kfo <- setdiff(getNames(cropAreaHr, dim = 1), kbe)
    cropTypesHr <- mbind(
      setNames(dimSums(cropAreaHr[, , kfo][, , "rainfed"], dim = 3), "crop.kfo.rainfed"),
      setNames(dimSums(cropAreaHr[, , kfo][, , "irrigated"], dim = 3), "crop.kfo.irrigated"),
      setNames(dimSums(cropAreaHr[, , kbe][, , "rainfed"], dim = 3), "crop.kbe.rainfed"),
      setNames(dimSums(cropAreaHr[, , kbe][, , "irrigated"], dim = 3), "crop.kbe.irrigated"),
      setNames(cropSplitHr[, , "crop_fallow"], "crop.fallow.total"),
      setNames(cropSplitHr[, , "crop_treecover"], "crop.treecover.total")
    )
    nonCropHr <- landHr[, , "crop", invert = TRUE]
    getNames(nonCropHr) <- paste0(getNames(nonCropHr), ".total")
    x <- mbind(cropTypesHr, nonCropHr)
    getSets(x, fulldim = FALSE) <- c("x.y.iso", "year", "land.type.water")
    toolExpectLessDiff(dimSums(x, dim = 3), dimSums(landHr, dim = 3), 1e-6,
                       "The land split covers all land of the grid cells")
    return(list(x = x, isocountries = FALSE, unit = "Mha", min = 0,
                description = "MAgPIE land by grid cell, land type, sub-type and water supply"))
  } else {
    stop("Can only report for outputFormat = ESM/ScenarioMIP/downscaledmagpie/downscaledmagpieSplit/",
         "downscaledmagpieCroparea")
  }
}
