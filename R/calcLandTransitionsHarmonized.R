#' calcLandTransitionsHarmonized
#'
#' Harmonisation of the land transitions of an input dataset (\code{\link{calcLandTransitionsInput}}) for
#' the downscaling along its transitions, the counterpart of \code{\link{calcLandHarmonized}} on the
#' transitions route.
#'
#' @details
#' With \code{harmonization = "none"} the transitions, states and start map of the input are passed on
#' unchanged: the downscaling starts from the start map of the input, and \code{target} and
#' \code{harmonizationPeriod} are not used. Other harmonisations are not available on the transitions
#' route yet. The designed approach harmonises the unit states of the main land types with
#' \code{\link{toolGetHarmonizer}} against \code{\link{calcLandTargetLowRes}}, gives the sub-types the
#' change of their main type in proportion to their area, and adds the change of the harmonisation offset
#' between two time steps as a stage of its own before the optimisation (from
#' \code{\link{toolMinimalTransitions}}), so that the transitions of the input stay as they are; the start
#' map is then the target map at the end of history.
#'
#' @inheritParams calcLandHarmonized
#' @param harmonization name of harmonization method; only "none" on the transitions route
#' @return list as returned by \code{\link{calcLandTransitionsInput}}
#' @author Patrick v. Jeetze
#' @seealso \code{\link{calcLandTransitionsHighRes}}
calcLandTransitionsHarmonized <- function(input, target, harmonizationPeriod, harmonization) {
  if (harmonization != "none") {
    stop("harmonization \"", harmonization, "\" is not available for downscaling = \"transitions\" yet, ",
         "only \"none\" (see the details of calcLandTransitionsHarmonized)")
  }
  x <- calcOutput("LandTransitionsInput", input = input, aggregate = FALSE)
  return(list(x = x, class = "list", unit = "Mha", cache = FALSE,
              description = "Land transitions of the input, harmonised"))
}
