#' calcBIIHighRes
#'
#' Biodiversity Intactness Index (BII) of each grid cell: the BII coefficients of the units weighted by the
#' shares of the BII land classes in the grid cell, from the land of \code{\link{calcLandTransitionsHighRes}}.
#'
#' @details input = "magpie": the BII coefficients of the clusters (\code{magpie4::BII}) by land class,
#' each forested and non-forested. Pasture is split into managed pasture and rangeland, and each land class
#' into forested and non-forested land, with the side layers (\code{luh2_side_layers_0.5.mz}); other land
#' into primary and secondary other land with the shares of the land-use initialisation
#' (\code{magpie4::PrimSecdOtherLand}). Primary other land has BII = 1.
#'
#' @inheritParams calcLandTransitionsHighRes
#' @return BII of each grid cell (unitless)
#' @author Patrick v. Jeetze
#' @seealso \code{\link{fullDOWNSCALEDMAGPIE}}
calcBIIHighRes <- function(input, harmonization = "none") {
  "!# @monitor magpie4:::PrimSecdOtherLand"

  landHr <- dimSums(calcOutput("LandTransitionsHighRes", input = input, harmonization = harmonization,
                               aggregate = FALSE)[["landHr"]], dim = 3.2)

  # BII coefficients (grid cells x classes) and the shares of the classes in the grid cells
  if (input == "magpie") {
    mapping <- readSource("MagpieFulldataGdx", subtype = "clustermap")
    biiLr <- readSource("MagpieFulldataGdx", subtype = "bii")
    sideLayersHr <- readSource("MagpieRunFiles", subtype = "sideLayers")
    landStartHr <- readSource("MagpieRunFiles", subtype = "landuseInit")[, "y1995", ]

    biiLr <- mbind(biiLr[, , "other", invert = TRUE],
                   setNames(biiLr[, , "other"], c("primother.forested", "primother.nonforested")),
                   setNames(biiLr[, , "other"], c("secdother.forested", "secdother.nonforested")))
    biiLr[, , c("primother.forested", "primother.nonforested")] <- 1
    coefficientHr <- toolAggregate(biiLr, mapping, from = "cluster", to = "cell")

    shareHr <- mbind(landHr[, , "past", invert = TRUE],
                     collapseNames(landHr[, , "past"]) * sideLayersHr[, , c("manpast", "rangeland")])
    shareHr <- shareHr / dimSums(shareHr, dim = 3)
    shareHr[is.nan(shareHr)] <- 0
    shareHr <- magpie4::PrimSecdOtherLand(shareHr, landStartHr)
    shareHr <- shareHr * sideLayersHr[, , c("forested", "nonforested")]
  } else {
    stop("Unsupported input dataset \"", input, "\"")
  }

  return(list(x = dimSums(shareHr * coefficientHr, dim = 3, na.rm = TRUE), isocountries = FALSE, unit = "1",
              description = "BII of each grid cell"))
}
