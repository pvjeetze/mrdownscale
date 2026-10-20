#' fullDOWNSCALEDMAGPIE
#'
#' Run the pipeline to generate downscaled MAgPIE data on the 0.5 degree grid of the land-use
#' initialisation (landuseinit). Full report on consistency checks is printed and written to
#' consistencyCheck.log.
#'
#' @details With \code{downscaling = "proportional"}, the land of the run is harmonised with
#' landuseinit and downscaled with \code{luscale::interpolate2}; the output is
#' \code{cell.land_0.5.mz} in the format of avl_land_t_0.5.mz.
#'
#' With \code{downscaling = "transitions"} (the default), the transitions of the run are allocated to
#' the grid cells (\code{\link{calcLandTransitionsHighRes}}) without harmonisation (\code{harmonization =
#' "none"}, the model starts from landuseinit). The product is selected with \code{product}:
#' \itemize{
#'   \item land: \code{cell.land_0.5.mz}, \code{cell.land_split_0.5.mz} (and their shares),
#'   \code{cell.croparea_0.5_share.mz}, \code{cell.land_transitions_0.5.mz}, \code{cell.peatland_0.5.mz}
#'   (and its share; from \code{s58_fix_peatland}), \code{cell.land_diagnostics.csv},
#'   \code{cell.land_constraint_excess.csv}, \code{cell.land_transition_change.csv},
#'   \code{cell.peatland_excess.csv};
#'   \item BII: \code{cell.bii_0.5.mz};
#'   \item landConservation: \code{cell.conservation_land_0.5.mz};
#'   \item ageClasses: \code{cell.ageclass_<type>_0.5_share.mz} for secdforest, forestry and other land, and
#'   \code{cell.ageclass_diagnostics.csv}.
#' }
#'
#' @param rev revision number of the data. If not provided the current date will be used instead.
#' When called via madrat::retrieveData rev will be converted to numeric_version.
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets (not used with harmonization "none")
#' @param target Name of dataset to be used as harmonization target and downscaling reference
#' @param downscaling name of downscaling method: "proportional" or "transitions", see details
#' @param harmonization name of harmonization method, see \code{\link{toolGetHarmonizer}}; "none" (only
#' with "transitions", the only method available there)
#' @param product the product to write with "transitions": "land", "BII", "landConservation" or
#' "ageClasses", see details
#'
#' @author Pascal Sauer, Patrick v. Jeetze
fullDOWNSCALEDMAGPIE <- function(rev = numeric_version("0"), # nolint: cyclocomp_linter
                                 harmonizationPeriod = c(2015, 2050), target = "landuseinit",
                                 downscaling = "transitions", harmonization = "none", product = "land") {
  if (target != "landuseinit") {
    stop("currently only target = landuseinit is supported")
  }

  # Writes a 0.5 degree output with the spatial sets of the outputs.
  .write <- function(x, file, comment) {
    message("Write outputs ", file)
    getSets(x, fulldim = FALSE)[1] <- "x.y.iso"
    write.magpie(x, file, comment = comment)
  }

  # Land report of the run in one of the MAgPIE formats.
  .report <- function(outputFormat) {
    return(calcOutput("LandReport", outputFormat = outputFormat, input = "magpie",
                      harmonizationPeriod = harmonizationPeriod, yearsSubset = seq(1995, 2100, 5),
                      harmonization = harmonization, downscaling = downscaling, aggregate = FALSE))
  }

  if (downscaling == "proportional") {
    if (product != "land") {
      stop("downscaling = \"proportional\" gives only the product \"land\"")
    }
    if (harmonization == "none") {
      stop("downscaling = \"proportional\" needs a harmonization method, see toolGetHarmonizer")
    }
    # LandReport with outputFormat = downscaledmagpie assumes target = landuseinit
    calcOutput("LandReport", outputFormat = "downscaledmagpie", input = "magpie",
               harmonizationPeriod = harmonizationPeriod, yearsSubset = seq(1995, 2100, 5),
               harmonization = harmonization, downscaling = downscaling,
               aggregate = FALSE, file = "cell.land_0.5.mz")
  } else if (downscaling != "transitions") {
    stop("Unsupported downscaling method \"", downscaling, "\"")
  } else if (harmonization != "none") {
    stop("downscaling = \"transitions\" is available only with harmonization = \"none\"")
  } else if (product == "land") {
    land <- calcOutput("LandTransitionsHighRes", input = "magpie", aggregate = FALSE)
    .write(land[["transitionsHr"]], "cell.land_transitions_0.5.mz",
           comment = "unit: Mha per grid-cell and time step, gross transitions by stage of the step")
    utils::write.csv(land[["diagnostics"]], "cell.land_diagnostics.csv", row.names = FALSE)
    if (!is.null(land[["excess"]])) {
      utils::write.csv(land[["excess"]], "cell.land_constraint_excess.csv", row.names = FALSE)
    }
    if (!is.null(land[["transitionChange"]])) {
      utils::write.csv(land[["transitionChange"]], "cell.land_transition_change.csv", row.names = FALSE)
    }

    if (!is.null(land[["peatlandHr"]])) {
      # peatland from the year in which the model starts to change it
      firstYear <- readSource("MagpieRunFiles", subtype = "config")[["gms"]][["s58_fix_peatland"]]
      peatHr <- land[["peatlandHr"]]
      peatHr <- peatHr[, getYears(peatHr, as.integer = TRUE) >= firstYear, ]
      utils::write.csv(land[["peatlandExcess"]], "cell.peatland_excess.csv", row.names = FALSE)
      .write(peatHr, "cell.peatland_0.5.mz", comment = "unit: Mha per grid-cell")
      # share of the grid cell area (area of a 0.5 degree grid cell in Mha)
      coords <- toolGetMappingCoord2Country(pretty = TRUE)
      cellAreaHr <- new.magpie(cells_and_regions = coords[["coords"]],
                               fill = 1e-10 * (111.263 * 1000 * 0.5)^2 * cos(coords[["lat"]] * pi / 180))
      peatShareHr <- peatHr / cellAreaHr
      peatShareHr[is.nan(peatShareHr) | is.infinite(peatShareHr)] <- 0
      .write(peatShareHr, "cell.peatland_0.5_share.mz", comment = "unit: grid-cell area fraction")
    }
    rm(land)

    .write(.report("downscaledmagpieCroparea"), "cell.croparea_0.5_share.mz",
           comment = "unit: croparea fractions of total grid-cell")
    landSplitHr <- .report("downscaledmagpieSplit")
    landHr <- .report("downscaledmagpie")
    .write(landSplitHr, "cell.land_split_0.5.mz", comment = "unit: Mha per grid-cell")
    .write(landSplitHr / dimSums(landSplitHr, dim = 3), "cell.land_split_0.5_share.mz",
           comment = "unit: grid-cell land area fraction")
    .write(landHr, "cell.land_0.5.mz", comment = "unit: Mha per grid-cell")
    .write(landHr / dimSums(landHr, dim = 3), "cell.land_0.5_share.mz", comment = "unit: grid-cell land area fraction")
  } else if (product == "BII") {
    .write(calcOutput("BIIHighRes", input = "magpie", aggregate = FALSE), "cell.bii_0.5.mz", comment = "unitless")
  } else if (product == "landConservation") {
    if (!"wdpa" %in% readSource("MagpieRunFiles", subtype = "available")) {
      stop("the run has no land conservation data (wdpa_baseline_0.5.mz)")
    }
    .write(calcOutput("LandConservationHighRes", input = "magpie", onLand = TRUE, aggregate = FALSE),
           "cell.conservation_land_0.5.mz", comment = "unit: Mha per grid-cell")
  } else if (product == "ageClasses") {
    ageClasses <- calcOutput("AgeClassesReport", outputFormat = "downscaledmagpie", input = "magpie",
                             aggregate = FALSE)
    diagnostics <- ageClasses[["diagnostics"]]
    utils::write.csv(diagnostics, "cell.ageclass_diagnostics.csv", row.names = FALSE)
    message("Largest deviation of the grid cell area from the grid cell land: ",
            signif(max(diagnostics[["gridCellDeviation"]]), 3),
            " Mha; summed deviation of the cluster age classes from the model: ",
            signif(max(diagnostics[["clusterDeviation"]]), 3), " Mha in a year and land type")
    # one file per land type: the memory use of write.magpie() grows with the number of items
    longName <- c(secdforest = "secondary forest", forestry = "forestry", other = "other land")
    for (type in names(longName)) {
      .write(ageClasses[[type]], paste0("cell.ageclass_", type, "_0.5_share.mz"),
             comment = paste0("unit: share of ", longName[[type]], " area in each age class"))
      ageClasses[[type]] <- NULL
      gc(verbose = FALSE)
    }
  } else {
    stop("Unsupported product \"", product, "\"")
  }

  toolWriteMadratLog(logPath = "consistencyCheck.log")
}
