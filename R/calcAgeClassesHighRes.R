#' calcAgeClassesHighRes
#'
#' Age classes of an input dataset on the grid, along the grid cell land and transitions of
#' \code{\link{calcLandTransitionsHighRes}}. In the first year each grid cell gets the age composition of its
#' unit; then the operations of \code{\link{calcAgeClassesInput}} are replayed in each time step: unit
#' cohort flows are allocated to the grid cells with \code{\link{toolAllocateCohortFlows}} (within the grid
#' cell totals of the land transitions where given), ageing with \code{\link{toolShiftCohorts}}. Area that
#' enters and leaves a land type in the same time step passes through without age. At the end of each time
#' step the grid cell stocks are fitted to the unit stocks with \code{\link{toolFitCohortStocks}}.
#'
#' @inheritParams calcLandTransitionsHighRes
#' @return List with magclass objects of the grid cell area by age class (Mha) for each group of
#' \code{\link{calcAgeClassesInput}}, and \code{diagnostics} (per year and land type: largest deviation of
#' the grid cell area from the grid cell land, summed deviation from the unit stocks, unit area, area passing
#' through, area moved by the fit)
#' @author Patrick v. Jeetze
#' @seealso \code{\link{calcAgeClassesReport}}
calcAgeClassesHighRes <- function(input, harmonization = "none") { # nolint: cyclocomp_linter
  # Unit values by age class in a year, without year.
  .unitAges <- function(x, year) {
    return(setYears(x[units, year, ageClasses], NULL))
  }

  # Adds amount (one value per grid cell) to the first nShift age classes of stock, in equal parts.
  .toEntry <- function(stock, amount, nShift) {
    entry <- ageClasses[seq_len(nShift)]
    stock[, , entry] <- stock[, , entry] + amount / nShift
    return(stock)
  }

  # Grid cell sum of the selected grid cell transitions of a year.
  .cellFlow <- function(flows, select) {
    if (!any(select)) {
      return(zeroHr)
    }
    return(dimSums(flows[, , select], dim = 3))
  }

  ages <- calcOutput("AgeClassesInput", input = input, aggregate = FALSE)
  land <- calcOutput("LandTransitionsHighRes", input = input, harmonization = harmonization, aggregate = FALSE)
  stocksLr <- ages[["stocksLr"]]
  pools <- ages[["pools"]]
  mapping <- ages[["mapping"]]
  years <- getYears(stocksLr)
  landHr <- land[["landHr"]][, years, ]
  transitionsHr <- land[["transitionsHr"]][, years, ]
  rm(land)

  units <- getItems(stocksLr, dim = 1)
  cells <- getCells(landHr)
  stopifnot("mapping lacks grid cells of the land" = all(cells %in% mapping[["cell"]]))
  ageClasses <- getItems(stocksLr, dim = "ac")
  yearNumbers <- getYears(stocksLr, as.integer = TRUE)
  zeroHr <- new.magpie(cells, NULL, NULL, fill = 0)

  # unit stock of each land type by age class, and its grid cell land
  stockLr <- list()
  landOfPool <- list()
  for (pool in pools) {
    stockLr[[pool]] <- collapseDim(stocksLr[, , pool], dim = c(3.1, 3.2))
    landOfPool[[pool]] <- collapseDim(landHr[, , pool], dim = 3)
  }

  # labels of the grid cell transitions
  stage <- getItems(transitionsHr, dim = 3.1, full = TRUE)
  from <- paste(getItems(transitionsHr, dim = 3.2, full = TRUE), getItems(transitionsHr, dim = 3.3, full = TRUE),
                sep = ".")
  to <- paste(getItems(transitionsHr, dim = 3.4, full = TRUE), getItems(transitionsHr, dim = 3.5, full = TRUE),
              sep = ".")

  # first year: the age composition of the unit in each grid cell (equal shares where the unit has no stock)
  stock <- list()
  for (pool in pools) {
    unitStock <- .unitAges(stockLr[[pool]], years[1])
    stockSum <- dimSums(unitStock, dim = 3)
    share <- unitStock / stockSum
    share[is.nan(share)] <- 1 / length(ageClasses)
    shareHr <- toolAggregate(share, mapping, from = "cluster", to = "cell")[cells, , ]
    stock[[pool]] <- shareHr * setYears(landOfPool[[pool]][, years[1], ], NULL)
  }

  # stocks of each group by year, joined once at the end: each assignment into a magclass object copies the
  # whole object, so filling one object year by year takes several times its size in memory
  out <- list()
  diagnostics <- list()

  for (i in seq_along(years)) {
    year <- years[i]
    passThrough <- stats::setNames(rep(0, length(pools)), pools)

    if (i > 1) {
      nShift <- (yearNumbers[i] - yearNumbers[i - 1]) / ages[["classWidth"]]
      stopifnot("time steps must be multiples of the age class width" = nShift == round(nShift))
      # grid cell transitions of the year
      flows <- setYears(transitionsHr[, year, ], NULL)

      for (operation in ages[["operations"]]) {
        if (operation[["kind"]] == "within") {
          pool <- operation[["pool"]]
          moved <- toolAllocateCohortFlows(stock[[pool]], .unitAges(operation[["flowLr"]], year), mapping)
          entry <- dimSums(moved, dim = 3)
          if (!is.null(operation[["inflowStages"]])) {
            entry <- entry + .cellFlow(flows, stage %in% operation[["inflowStages"]] & to == pool)
          }
          stock[[pool]] <- .toEntry(stock[[pool]] - moved, entry, nShift)
        } else if (operation[["kind"]] == "transfer") {
          cellTotal <- .cellFlow(flows, stage %in% operation[["stages"]] & from == operation[["from"]] &
                                   to == operation[["to"]])
          moved <- toolAllocateCohortFlows(stock[[operation[["from"]]]], .unitAges(operation[["flowLr"]], year),
                                           mapping, cellTotal, tol = 1e-6)
          stock[[operation[["from"]]]] <- stock[[operation[["from"]]]] - moved
          stock[[operation[["to"]]]] <- stock[[operation[["to"]]]] + moved
        } else if (operation[["kind"]] == "shift") {
          stock <- lapply(stock, toolShiftCohorts, nShift = nShift)
        } else if (operation[["kind"]] == "exchange") {
          pool <- operation[["pool"]]
          outside <- !stage %in% operation[["excludeStages"]]
          outflow <- .cellFlow(flows, outside & from == pool & to != pool)
          inflow <- .cellFlow(flows, outside & to == pool & from != pool)
          # area that enters and leaves in the same time step passes through without age
          through <- pmin(pmax(outflow - dimSums(stock[[pool]], dim = 3), 0), inflow)
          passThrough[[pool]] <- sum(through)

          left <- toolAllocateCohortFlows(stock[[pool]], .unitAges(operation[["outflowLr"]], year), mapping,
                                          outflow - through, tol = 1e-6)
          stock[[pool]] <- stock[[pool]] - left
          regrowth <- toolAllocateCohortFlows(stock[[pool]], .unitAges(operation[["regrowthLr"]], year), mapping)
          stock[[pool]] <- .toEntry(stock[[pool]] - regrowth, dimSums(regrowth, dim = 3) + inflow - through, nShift)
        } else {
          stop("Unsupported operation \"", operation[["kind"]], "\"")
        }
      }
    }

    # fit to the unit stocks, store and diagnose
    for (pool in pools) {
      unitStock <- .unitAges(stockLr[[pool]], year)
      tracked <- stock[[pool]]
      stock[[pool]] <- toolFitCohortStocks(tracked, unitStock, mapping, tol = 1e-6)
      type <- sub("^.*[.]", "", pool)
      fittedLr <- toolAggregate(stock[[pool]], mapping, from = "cell", to = "cluster")[units, , ]
      diagnostics[[length(diagnostics) + 1]] <- data.frame(
        year = year, type = ifelse(type == "total", sub("[.].*$", "", pool), type),
        gridCellDeviation = max(abs(dimSums(stock[[pool]], dim = 3) - setYears(landOfPool[[pool]][, year, ], NULL))),
        clusterDeviation = sum(abs(fittedLr - unitStock)),
        clusterArea = sum(unitStock), passThrough = passThrough[[pool]],
        movedByFit = sum(abs(stock[[pool]] - tracked)) / 2
      )
    }
    for (group in names(ages[["groups"]])) {
      total <- 0
      for (pool in ages[["groups"]][[group]]) {
        total <- total + stock[[pool]]
      }
      out[[group]][[year]] <- setYears(total, year)
    }
  }

  for (group in names(out)) {
    out[[group]] <- mbind(out[[group]])
    getSets(out[[group]], fulldim = FALSE) <- c("x.y.iso", "year", "ac")
  }
  out[["diagnostics"]] <- do.call(rbind, diagnostics)
  return(list(x = out, class = "list", unit = "Mha",
              description = "Age classes of the input by grid cell"))
}
