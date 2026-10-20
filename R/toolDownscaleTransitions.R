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
  units <- getItems(x, 1)
  years <- getYears(x)[-1]

  # unit states at the end of each time step (units x types x years)
  statesLr <- array(0, c(length(units), length(types), length(years)), dimnames = list(units, types, years))
  for (t in years) {
    statesLr[, , t] <- matrix(x[, t, types], nrow = length(units))
  }

  # net transitions of each time step: the annual transitions times the length of the time step,
  # labelled with the end year
  annual <- toolMinimalTransitions(x[, , types])
  stepLength <- diff(getYears(x, as.integer = TRUE))
  transitionsLr <- list()
  for (i in seq_along(years)) {
    total <- matrix(annual[, i, ], nrow = length(units), dimnames = list(units, getItems(annual, 3))) * stepLength[i]
    transitionsLr[[years[i]]] <- list(change = total[, colSums(total) > 0, drop = FALSE])
  }

  startHr <- matrix(xTarget[, , types], nrow = ncells(xTarget), dimnames = list(getItems(xTarget, 1), types))
  result <- toolAllocateTransitions(startHr, mapping, statesLr, transitionsLr)

  out <- new.magpie(rownames(startHr), years, types, fill = 0)
  for (t in years) {
    out[, t, ] <- result[["stateHr"]][, , t]
  }
  getSets(out) <- c("x", "y", "year", "data")
  return(out)
}
