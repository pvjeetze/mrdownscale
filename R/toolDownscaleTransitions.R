#' toolDownscaleTransitions
#'
#' Downscaling along transitions: the net transitions between the states of each low resolution unit
#' (\code{\link{toolMinimalTransitions}}) are allocated to its grid cells with
#' \code{\link{toolAllocateTransitions}}, starting from the target map. The grid cell states sum to the
#' unit states in every year, and the grid cell transitions to the unit transitions.
#'
#' @param x magclass containing land to be downscaled, the first year being the year of \code{xTarget}
#' @param xTarget magclass target land use dataset in the first year of \code{x}
#' @param mapping mapping between \code{x} and \code{xTarget} (columns cell and lowRes)
#' @return downscaled land use dataset in the years of \code{x} after the first
#' @author Patrick v. Jeetze
#' @seealso \code{\link{toolDownscaleSimpleProportional}}
toolDownscaleTransitions <- function(x, xTarget, mapping) {
  mapping <- data.frame(cell = mapping[["cell"]], cluster = mapping[["lowRes"]])
  stopifnot(setequal(mapping[["cluster"]], getItems(x, 1)),
            setequal(mapping[["cell"]], getItems(xTarget, 1)),
            nyears(xTarget) == 1,
            nyears(x) >= 2,
            setequal(getItems(x, 3), getItems(xTarget, 3)))
  types <- getItems(x, 3)
  years <- getYears(x)[-1]

  # net transitions of each time step: the annual transitions times the length of the time step,
  # labelled with the end year, as one stage
  annual <- toolMinimalTransitions(x[, , types])
  stepLength <- diff(getYears(x, as.integer = TRUE))
  transitionsLr <- annual
  for (i in seq_along(years)) {
    transitionsLr[, i, ] <- annual[, i, ] * stepLength[i]
  }
  transitionsLr <- add_dimension(setYears(transitionsLr, years), dim = 3.1, add = "stage", nm = "change")

  result <- toolAllocateTransitions(xTarget[, , types], mapping, x[, years, types], transitionsLr)
  out <- result[["stateHr"]]
  getSets(out) <- c("x", "y", "year", "data")
  return(out)
}
