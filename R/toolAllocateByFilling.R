#' Allocate by filling
#'
#' Distributes totals of coarse units (Lr, for example the clusters of a model) to their grid cells (Hr) in
#' proportion to a prior, without exceeding the cap of any grid cell. Where the grid cells with a prior
#' are full, the rest goes to the other grid cells of the unit in proportion to their remaining cap.
#'
#' @details For each unit, year and item, the grid cell values are \code{min(cap, t * prior)} with the
#' factor \code{t} chosen so that they sum to the unit total (closed form, no iteration). Where the
#' caps of the grid cells with a prior are too small, these grid cells get their cap and the rest is
#' distributed in the same way over the remaining cap of all grid cells of the unit. Where the caps
#' of all grid cells are too small, every grid cell gets its cap, and the unit total is not reached.
#' As long as the grid cells with a prior can hold the total, this is the solution of
#' \code{\link{toolAllocateByCrossEntropy}} with one hard cap per grid cell and item, computed directly. The
#' spill-over to grid cells without a prior has no counterpart in \code{\link{toolAllocateByCrossEntropy}}, where a
#' zero prior stays zero; use \code{\link{toolAllocateByCrossEntropy}} for caps shared by several flows, soft or
#' signed caps.
#'
#' @param totalLr magclass object of unit totals (units x years x items).
#' @param priorHr magclass object of prior weights (grid cells x the years and items of \code{totalLr}).
#' @param capHr magclass object of caps with the dimensions and names of \code{priorHr}. Grid cells with a
#' prior or a cap of zero or below (for example -1e-17 from rounding) take nothing.
#' @inheritParams toolAllocateByCrossEntropy
#' @return magclass object of the grid cell values with the dimensions and names of \code{priorHr}.
#' @author Patrick v. Jeetze
#' @seealso \code{\link{toolAllocateByCrossEntropy}}
#' @examples
#' total <- new.magpie("A", NULL, "protected", fill = 1.5)
#' prior <- new.magpie(c("c1", "c2", "c3"), NULL, "protected", fill = c(0.2, 0.1, 0))
#' cap <- new.magpie(c("c1", "c2", "c3"), NULL, "protected", fill = 1)
#' mapping <- data.frame(cell = c("c1", "c2", "c3"), cluster = "A")
#' # c1 and c2 are filled first (c1 reaches its cap), the rest goes to c3
#' toolAllocateByFilling(total, prior, cap, mapping)
#' @export
toolAllocateByFilling <- function(totalLr, priorHr, capHr, mapping) {
  hrToLr <- as.character(mapping[["cluster"]][match(getCells(priorHr), mapping[["cell"]])])
  stopifnot("totalLr, priorHr and capHr must be magclass objects" =
              is.magpie(totalLr) & is.magpie(priorHr) & is.magpie(capHr),
            "priorHr needs the grid cell names of mapping" = !anyNA(hrToLr),
            "totalLr needs all units of the grid cells" = all(hrToLr %in% getCells(totalLr)),
            "priorHr needs the years and items of totalLr" =
              identical(getYears(priorHr), getYears(totalLr)) & identical(getNames(priorHr), getNames(totalLr)),
            "capHr needs the dimensions and names of priorHr" = identical(dimnames(capHr), dimnames(priorHr)))

  # Grid cell values of one unit: min(cap, t * prior) with sum total, or the cap where the caps are too
  # small.
  .fillCells <- function(prior, cap, total) {
    x <- rep(0, length(prior))
    free <- prior > 0 & cap > 0
    remaining <- total
    while (any(free) && remaining > 0) {
      proposal <- remaining / sum(prior[free]) * prior
      full <- free & proposal >= cap
      if (!any(full)) {
        x[free] <- proposal[free]
        break
      }
      x[full] <- cap[full]
      remaining <- remaining - sum(cap[full])
      free <- free & !full
    }
    return(x)
  }

  # the units are filled one by one on arrays: assigning the values of each unit to a magclass object in this
  # loop takes about 4000 times as long
  prior <- as.array(priorHr)
  cap <- as.array(capHr)
  total <- as.array(totalLr)
  out <- prior * 0
  cellsOf <- split(seq_along(hrToLr), hrToLr)
  for (t in seq_len(nyears(totalLr))) {
    for (item in seq_len(ndata(totalLr))) {
      for (unit in names(cellsOf)) {
        cells <- cellsOf[[unit]]
        x <- .fillCells(prior[cells, t, item], cap[cells, t, item], total[unit, t, item])
        rest <- total[unit, t, item] - sum(x)
        if (rest > 0) {
          room <- cap[cells, t, item] - x
          x <- x + .fillCells(room, room, rest)
        }
        out[cells, t, item] <- x
      }
    }
  }

  result <- priorHr
  result[] <- out
  return(result)
}
