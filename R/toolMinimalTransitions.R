#' toolMinimalTransitions
#'
#' Transitions between land categories from the change of their states (a land data set with at least
#' 2 time steps), for data with states only. The net loss of each category goes to the categories with
#' a net gain, in proportion to their gain; this is the smallest gross change that gives the states
#' (equal preference of transitions across all categories). Transitions to primf and primn are removed.
#'
#' In addition to the net effect it can also estimate gross transition. For that
#' purpose a reference data set containing bidirectional transition shares
#' must be provided.
#'
#' If the time step length is longer than 1 year the returned object contains
#' reference years for each period which can be repeated to retrieve the
#' full transition between two time periods, e.g. if you provide two time steps
#' 2000 and 2005 the return value will be the transition for year 2001. Repeating
#' the same transition also in 2002, 2003, 2004 and 2005 will give the full
#' transition from 2000 to 2005.
#'
#' @param x magpie data set containing land data with at least two time steps
#' to extract net transitions from
#' @param gross either boolean or a magpie object containing bidirectional
#' transition shares relative to the area of the involved land pools (transition
#' divided by the area of the land pool in the "from" sub dimension). If set to
#' FALSE only net transitions will be returned. If set to TRUE an internal
#' gross transition estimate based on average gross transitions in LUH2 in the
#' period from 1995 to 2015 will be used.
#' @author Jan Philipp Dietrich, Patrick v. Jeetze
#' @seealso \code{\link{toolDownscaleTransitions}}
#' @export

toolMinimalTransitions <- function(x, gross = FALSE) {
  stopifnot(nyears(x) >= 2)
  # assign transitions of a period to the first year in this transition period
  years <- paste0("y", getYears(x, as.integer = TRUE)[1:(nyears(x) - 1)] + 1)
  diff <- setItems(x[, 2:nyears(x), ], years, dim = 2) - setItems(x[, 1:(nyears(x) - 1), ], years, dim = 2)
  reduce <- -diff
  reduce[reduce < 0] <- 0 # reduce is all positive values, so needs to be subtracted, not added
  expand <- diff
  expand[expand < 0] <- 0
  stopifnot(identical(diff, expand - reduce))

  split <- expand / dimSums(expand, dim = 3) # share of total expansion per landtype
  split[is.na(split)] <- 0
  stopifnot(all.equal(diff, split * dimSums(expand, dim = 3) - reduce, check.attributes = FALSE))

  rm(expand, diff)
  gc()

  getSets(reduce, fulldim = FALSE)[3] <- "from"
  getSets(split, fulldim = FALSE)[3] <- "to"

  withr::with_options(list(magclass_setMatching = TRUE, magclass_sizeLimit = 10e+10), {
    out <- reduce * split
  })

  # remove entries for transitions to itself (1) as well as transitions to primn
  # or primf (2) as they are not allowed and should not happen
  primfn <- grep("\\.(primn|primf)$", getItems(out, dim = 3), value = TRUE)
  if (any(out[, , primfn] != 0)) {
    toolStatusMessage("warn", paste("removing non-zero transitions to primf/primn,",
                                    "transitions are likely inconsistent to states"))
  }
  transitionsToRemove <- union(paste0(getItems(reduce, dim = 3), ".", getItems(reduce, dim = 3)), primfn)
  out <- out[, , transitionsToRemove, invert = TRUE]

  # correct for timestep length (transitions should reflect a single representative
  # year, not the full transition.)
  tLengths <- new.magpie(years = getYears(out))
  tLengths[, , ] <- diff(getYears(x, as.integer = TRUE))
  if (any(tLengths != 1)) out <- out / tLengths

  if (!isFALSE(gross)) {
    if (isTRUE(gross)) {
      gross <- read.magpie(system.file("extdata/meanBidirectionalTransitionsShares1995to2015.mz",
                                       package = "mrdownscale"))
    }
    if (!is.magpie(gross)) stop('"gross" must be a MAgPIE object containing bidirectional transition shares!')

    # add (empty) self-transitions to enable following computations
    selfTransitions <- paste0(getItems(gross, dim = 3.1), ".", getItems(gross, dim = 3.1))
    tmp <- gross[, , 1]
    tmp[, , ] <- 0
    tmp <- setItems(tmp[, , rep(1, length(selfTransitions))], selfTransitions, dim = 3, raw = TRUE)
    gross <- mbind(gross, tmp)
    rm(tmp)
    getSets(x)[4] <- "from"
    withr::with_options(list(magclass_setMatching = TRUE, magclass_sizeLimit = 10e+10), {
      x2 <- magpie_expand(x[, , getItems(gross, dim = "from")], gross[getItems(x, dim = 1), , ])
      gross2 <- magpie_expand(gross[getItems(x, dim = 1), , ], x2)
      getSets(gross2)[3] <- "year"
      grossOptions <- x2 * gross2
    })
    rm(x2, gross, gross2)

    # remove irrelevant entries
    grossOptions <- grossOptions[, , selfTransitions, invert = TRUE]

    # grossOptions contains now 4 (2x2) options for gross transitions for every
    # connection: values computed as shares relative to each both involved landtypes
    # computed relative to the land at start and end of the transition
    # (2 land types x 2 time steps)
    # in the following always the smallest of these 4 values is being selected
    # as gross transition to prevent overbooking of land areas
    invItems <- sub("^(.*)\\.(.*)$", "\\2.\\1", getItems(grossOptions, dim = 3))
    years <- paste0("y", getYears(grossOptions, as.integer = TRUE)[1:(dim(grossOptions)[2] - 1)] + 1)
    opt1 <- setItems(grossOptions[, 2:dim(grossOptions)[2], ], dim = 2, years)
    opt2 <- setItems(grossOptions[, 1:(dim(grossOptions)[2] - 1), ], dim = 2, years)[, , getItems(opt1, dim = 3)]
    opt3 <- setItems(opt1, dim = 3, invItems, raw = TRUE)[, , getItems(opt1, dim = 3)]
    opt4 <- setItems(opt2, dim = 3, invItems, raw = TRUE)[, , getItems(opt1, dim = 3)]
    gross <- pmin(opt1, opt2, opt3, opt4)
    out[, , getItems(gross, dim = 3)] <- out[, , getItems(opt1, dim = 3)] + gross
  }

  return(out)
}
