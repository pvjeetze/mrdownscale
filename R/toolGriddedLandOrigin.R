#' Track the origin of gridded land
#'
#' Area of land types in each grid cell by the type it came from, along gross grid cell transitions (as from
#' \code{\link{calcLandTransitionsHighRes}}). Area that enters a tracked type from another type takes that type
#' as origin, area that moves between tracked types keeps its origin and its year of entry, and the initial
#' area of a type has the type itself as origin. Outflow from a tracked type takes the area that entered last
#' first (\code{"newest"}), the area that entered first first (\code{"oldest"}), or takes from all area in
#' proportion (\code{"proportional"}).
#'
#' @param transitionsHr magclass object of the gross grid cell transitions per year with names
#' "stage.from.to". The stages of a year are replayed in the order of their items. The outflow of a
#' stage comes from the area at the start of the stage, and only beyond it from the inflow of the
#' stage (area passing through).
#' @param startHr magclass object of the grid cell area by type in the year before the first year of
#' \code{transitionsHr} (one year, with the grid cells of \code{transitionsHr}).
#' @param types Types to track (items of "from" and "to" in \code{transitionsHr}).
#' @param order Which area an outflow takes: "newest", "oldest" or "proportional", see Description. One value
#' for all outflows, or a named vector by transition ("from.to") with an element "default" for the
#' transitions not named (a named transition need not occur in \code{transitionsHr}). Where the outflows
#' of a type in one stage take in different orders, the oldest area is taken first, then the newest, then
#' in proportion, and outflows of one order with different \code{last} in turn; each also from the inflow
#' of the stage passing through.
#' @param last Origins that an outflow takes only after all other area of the type in the grid cell (for
#' example protected land), in the order of \code{order}: a list by transition ("from.to") of character
#' vectors with an element "default" for the transitions not named, as \code{order}.
#' @return magclass object of the grid cell area by type and origin (names "type.origin") in the year of
#' \code{startHr} and the years of \code{transitionsHr}, with the attribute "takenLast": the area of the
#' origins in \code{last} taken by the outflows of each year, summed over the grid cells.
#' @author Patrick v. Jeetze
#' @seealso \code{\link{calcLandTransitionsHighRes}}
#' @examples
#' initial <- new.magpie("A.1", "y2000", c("crop", "other"), c(2, 1))
#' transitions <- new.magpie("A.1", c("y2010", "y2020"),
#'                           c("change.crop.other", "change.other.crop"), c(0.5, 0, 0, 0.7))
#' getSets(transitions, fulldim = FALSE)[3] <- "stage.from.to"
#' toolGriddedLandOrigin(transitions, initial, "other")
#' toolGriddedLandOrigin(transitions, initial, "other",
#'                       order = c(default = "newest", other.crop = "oldest"))
#' @export
toolGriddedLandOrigin <- function(transitionsHr, startHr, types, order = "newest", # nolint: cyclocomp_linter
                                  last = list(default = character(0))) {
  # The area by origin and year of entry is kept on arrays: the outflow takes the entries one by one, and
  # assigning these slices to a magclass object in a loop takes about 4000 times as long.

  # One stage of toolGriddedLandOrigin (area: list by tracked type; flows: grid cells x transitions of the stage;
  # orderOf, lastOf: the order and the origins taken last of each transition). The outflow comes from the area at
  # the start of the stage, and beyond it from the inflow of the stage (area passing through); it is carried
  # to tracked types with its origin, the area at the start with its year of entry. Inflow enters at entry.
  .originStage <- function(area, flows, from, to, orderOf, lastOf, origins, entry) {
    types <- names(area)
    nCell <- nrow(flows)

    # inflow from types that are not tracked, by origin (grid cells x origins)
    inflow <- list()
    arriving <- list()
    carried <- list()
    for (j in types) {
      incoming <- which(to == j & !from %in% types)
      inflow[[j]] <- flows[, incoming, drop = FALSE] %*% outer(from[incoming], origins, "==")
      arriving[[j]] <- matrix(0, nCell, length(origins))
      carried[[j]] <- 0
    }

    shortfall <- 0
    takenLast <- 0
    group <- paste(orderOf, vapply(lastOf, paste, character(1), collapse = ","))
    for (k in intersect(types, from)) {
      # the outflow of each order (and origins taken last) in turn: from the area at the start of the stage,
      # beyond it from the inflow of the stage
      for (o in intersect(c("oldest", "newest", "proportional"), orderOf[from == k])) {
        for (g in unique(group[from == k & orderOf == o])) {
          leaving <- which(from == k & group == g)
          total <- rowSums(flows[, leaving, drop = FALSE])
          taken <- .takeArea(area[[k]], total, o, which(origins %in% lastOf[[leaving[1]]]))
          area[[k]] <- taken[["area"]]
          takenLast <- takenLast + taken[["takenLast"]]
          inflowSum <- rowSums(inflow[[k]])
          passing <- inflow[[k]] * pmin(taken[["shortfall"]] / replace(inflowSum, inflowSum <= 0, 1), 1)
          inflow[[k]] <- inflow[[k]] - passing
          shortfall <- max(shortfall, taken[["shortfall"]] - rowSums(passing))
          for (m in leaving[to[leaving] %in% types]) {
            share <- flows[, m] / replace(total, total <= 0, 1)
            carried[[to[m]]] <- carried[[to[m]]] + taken[["taken"]] * share
            arriving[[to[m]]] <- arriving[[to[m]]] + passing * share
          }
        }
      }
    }

    for (j in types) {
      area[[j]] <- area[[j]] + carried[[j]]
      area[[j]][, , entry] <- area[[j]][, , entry] + inflow[[j]] + arriving[[j]]
    }
    return(list(area = area, shortfall = shortfall, takenLast = takenLast))
  }

  # Takes the grid cell outflow total from area (grid cells x origins x years of entry) as .takeInOrder, the
  # origins last (indices) only beyond all other area; returns also the area taken of these origins.
  .takeArea <- function(area, total, order, last) {
    if (length(last) == 0) {
      return(c(.takeInOrder(area, total, order), list(takenLast = 0)))
    }
    other <- area
    other[, last, ] <- 0
    first <- .takeInOrder(other, total, order)
    second <- .takeInOrder(area - other, first[["shortfall"]], order)
    taken <- first[["taken"]] + second[["taken"]]
    return(list(area = area - taken, taken = taken, shortfall = second[["shortfall"]],
                takenLast = sum(second[["taken"]])))
  }

  # Takes the grid cell outflow total from area (grid cells x origins x years of entry), the newest or the
  # oldest entries first or in proportion; returns the remaining area, the area taken and the outflow above
  # the area of each grid cell.
  .takeInOrder <- function(area, total, order) {
    stock <- rowSums(area)
    if (order == "proportional") {
      share <- matrix(pmin(total / replace(stock, stock <= 0, 1), 1), nrow(area), dim(area)[3])
    } else {
      # the share taken from the area of each year of entry (grid cells x years of entry), in turn
      entryArea <- rowSums(aperm(area, c(1, 3, 2)), dims = 2)
      share <- entryArea * 0
      remaining <- total
      for (e in sort(seq_len(ncol(entryArea)), decreasing = order == "newest")) {
        take <- pmin(remaining, entryArea[, e])
        share[, e] <- take / replace(entryArea[, e], entryArea[, e] <= 0, 1)
        remaining <- remaining - take
      }
    }
    taken <- area * as.vector(share[, rep(seq_len(ncol(share)), each = ncol(area))])
    return(list(area = area - taken, taken = taken, shortfall = pmax(total - stock, 0)))
  }

  # Area by grid cell and origin of each tracked type (summed over the years of entry), as one matrix
  # grid cells x (types and origins).
  .originArea <- function(area) {
    return(do.call(cbind, lapply(area, rowSums, dims = 2)))
  }

  # TRUE if the names of a setting by transition are "default" once and otherwise of the form "from.to"
  .byTransition <- function(setting) {
    return(sum(names(setting) == "default") == 1 &
             all(grepl("^[^.]+[.][^.]+$", setdiff(names(setting), "default"))))
  }

  stage <- getItems(transitionsHr, "stage", full = TRUE)
  from <- getItems(transitionsHr, "from", full = TRUE)
  to <- getItems(transitionsHr, "to", full = TRUE)
  if (is.null(names(order))) {
    order <- c(default = order)
  }
  stopifnot("order must be \"newest\", \"oldest\" or \"proportional\"" =
              all(order %in% c("newest", "oldest", "proportional")),
            "order needs one element \"default\", the other names are transitions \"from.to\"" =
              .byTransition(order),
            "last needs one element \"default\", the other names are transitions \"from.to\"" =
              is.list(last) & .byTransition(last),
            "startHr must have one year" = nyears(startHr) == 1,
            "grid cells of startHr and transitionsHr differ" = identical(getCells(startHr), getCells(transitionsHr)))

  transition <- paste(from, to, sep = ".")
  orderOf <- unname(order[transition])
  orderOf[is.na(orderOf)] <- order[["default"]]
  lastOf <- unname(last[transition])
  lastOf[!transition %in% names(last)] <- list(last[["default"]])
  years <- getYears(transitionsHr)
  origins <- unique(c(types, from[to %in% types & !from %in% types]))
  nCell <- ncells(startHr)

  # area by grid cell, origin and year of entry, for each tracked type
  area <- lapply(stats::setNames(types, types), function(k) {
    x <- array(0, c(nCell, length(origins), length(years) + 1))
    x[, match(k, origins), 1] <- as.vector(startHr[, , k])
    return(x)
  })

  out <- array(0, c(nCell, length(years) + 1, length(types) * length(origins)),
               dimnames = list(getCells(startHr), c(getYears(startHr), years),
                               paste(rep(types, each = length(origins)), origins, sep = ".")))
  out[, 1, ] <- .originArea(area)

  shortfall <- 0
  takenLast <- stats::setNames(rep(0, length(years)), years)
  for (i in seq_along(years)) {
    flowsYear <- matrix(transitionsHr[, years[i], ], nrow = nCell)
    for (s in unique(stage)) {
      n <- which(stage == s & from != to & (from %in% types | to %in% types))
      step <- .originStage(area, flowsYear[, n, drop = FALSE], from[n], to[n], orderOf[n], lastOf[n], origins,
                           i + 1)
      area <- step[["area"]]
      shortfall <- max(shortfall, step[["shortfall"]])
      takenLast[i] <- takenLast[i] + step[["takenLast"]]
    }
    out[, i + 1, ] <- .originArea(area)
  }

  if (shortfall > 1e-6) {
    warning("outflow exceeds the tracked area by up to ", signif(shortfall, 3), " in a grid cell")
  }

  result <- as.magpie(out, spatial = 1, temporal = 2)
  getSets(result, fulldim = FALSE) <- c(getSets(startHr, fulldim = FALSE)[1], "year", "type.origin")
  attr(result, "takenLast") <- takenLast
  return(result)
}
