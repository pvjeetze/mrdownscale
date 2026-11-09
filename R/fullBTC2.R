#' fullBTC2
#'
#' The land of the downscaling along the transitions in the land classes of the project Bending the Curve
#' 2.0 (\code{\link{calcLandReportBtC2}}), as grid cell products: \code{cell.land_btc2_0.5.mz} (Mha by grid
#' cell and class) and a NetCDF file with the share of each class in the grid cell area
#' (\code{LC_area_share}, classes coded 1 to 10 in the order of the classes) and the grid cell area
#' (\code{pixel_area}, Mha).
#'
#' @param rev revision number of the data. If not provided the current date will be used instead.
#' @param projectName prefix of the NetCDF file name (left out where the title of the run contains it)
#' @param logFile name of the log of the consistency checks
#' @author Patrick v. Jeetze
#' @seealso \code{\link{downscaleRun}}
fullBTC2 <- function(rev = numeric_version("0"), projectName = "BtC2v04", logFile = "consistencyCheck_BtC2.log") {
  landHr <- calcOutput("LandReportBtC2", input = "magpie", aggregate = FALSE)
  write.magpie(landHr, "cell.land_btc2_0.5.mz", comment = "unit: Mha per grid-cell")

  cfg <- readSource("MagpieRunFiles", subtype = "config")
  prefix <- if (grepl(tolower(projectName), tolower(cfg$title))) "" else paste0(projectName, "_")
  classes <- getItems(landHr, dim = 3)
  classCodes <- paste(seq_along(classes), classes, sep = "=")

  # grid cells on the 0.5 degree grid, the latitude from north to south
  lon <- seq(-179.75, 179.75, by = 0.5)
  lat <- -seq(-89.75, 89.75, by = 0.5)
  coords <- getCoords(landHr)
  index <- cbind(match(coords[["x"]], lon), match(coords[["y"]], lat))
  landArea <- dimSums(landHr, dim = 3)
  shareHr <- landHr / landArea

  shareArray <- array(NA, c(length(lon), length(lat), length(classes), nyears(landHr)))
  for (k in seq_along(classes)) {
    for (i in seq_len(nyears(landHr))) {
      shareArray[cbind(index, k, i)] <- as.vector(shareHr[, i, k])
    }
  }
  areaArray <- array(NA, c(length(lon), length(lat)))
  areaArray[index] <- as.vector(landArea[, 1, ])

  dimLon <- ncdf4::ncdim_def("lon", "degrees_east", lon)
  dimLat <- ncdf4::ncdim_def("lat", "degrees_north", lat)
  dimTime <- ncdf4::ncdim_def("time", "years", as.numeric(getYears(landHr, as.integer = TRUE)), calendar = "standard",
                              unlim = TRUE)
  dimClass <- ncdf4::ncdim_def("lc_class", paste(classCodes, collapse = "/"), seq_along(classes))
  shareVar <- ncdf4::ncvar_def("LC_area_share", units = "share of pixel area",
                               longname = "share of pixel occupied by various land covers",
                               dim = list(dimLon, dimLat, dimClass, dimTime), missval = NaN, prec = "double",
                               compression = 9)
  areaVar <- ncdf4::ncvar_def("pixel_area", units = "million ha", longname = "total area of the pixel",
                              dim = list(dimLon, dimLat), missval = NaN, prec = "double", compression = 9)
  nc <- ncdf4::nc_create(paste0(prefix, cfg$title, ".nc"), list(shareVar, areaVar), force_v4 = TRUE)
  ncdf4::ncvar_put(nc, areaVar, areaArray)
  ncdf4::ncvar_put(nc, shareVar, shareArray)
  ncdf4::ncatt_put(nc, 0, "title", "MAgPIE land cover projections")
  ncdf4::ncatt_put(nc, 0, "scenario", cfg$title)
  ncdf4::ncatt_put(nc, 0, "institution", "PIK")
  ncdf4::ncatt_put(nc, 0, "contact",
                   "F. Humpen\u00f6der, humpenoeder@pik-potsdam.de; P. v. Jeetze, vjeetze@pik-potsdam.de")
  ncdf4::ncatt_put(nc, 0, "history", date())
  ncdf4::nc_close(nc)

  toolWriteMadratLog(logPath = logFile)
}
