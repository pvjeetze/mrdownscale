#' toolDownscaleSimpleProportional
#'
#' Proportional downscaling with \code{luscale::interpolate2}: in each low resolution unit, a land type
#' that shrinks shrinks in each grid cell by the same share, and the freed area goes to the land types
#' that grow, in proportion to their growth.
#'
#' @param x magclass containing land to be downscaled
#' @param xTarget magclass target land use dataset for initialization year
#' @param xTargetLowRes magclass target land use dataset for initialization year in low resolution (like x)
#' @param mapping mapping between \code{x} and \code{xTarget}
#' @return downscaled land use dataset
#' @author Jan Philipp Dietrich, Pascal Sauer
toolDownscaleSimpleProportional <- function(x, xTarget, xTargetLowRes, mapping) {
  mapping$cluster <- mapping$lowRes
  mapping <- mapping[, c("cell", "cluster")]

  stopifnot(setequal(mapping$cluster, getItems(x, 1)),
            setequal(mapping$cell, getItems(xTarget, 1)),
            nyears(xTarget) == 1,
            getYears(xTargetLowRes) == getYears(xTarget))

  "!# @monitor luscale::interpolate2"

  # interpolate2 assumes constant total over time, but only warns if unfulfilled, convert that to error
  tryCatch({
    out <- luscale::interpolate2(x = x[, -1, ],
                                 x_ini = xTarget,
                                 x_ini_lr = xTargetLowRes,
                                 map = mapping)
  }, warning = stop)
  getSets(out)[1:2] <- c("x", "y")

  return(out)
}
