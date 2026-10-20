#' Allocate by filling
#'
#' Distributes totals of coarse units (Lr, for example the clusters of a model) to their grid cells (Hr) in
#' proportion to a prior, without exceeding the cap of any grid cell. Where the grid cells with a prior
#' are full, the rest goes to the other grid cells of the unit in proportion to their remaining cap.
#'
#' @details For each unit and column, the grid cell values are \code{min(cap, t * prior)} with the
#' factor \code{t} chosen so that they sum to the unit total (closed form, no iteration). Where the
#' caps of the grid cells with a prior are too small, these grid cells get their cap and the rest is
#' distributed in the same way over the remaining cap of all grid cells of the unit. Where the caps
#' of all grid cells are too small, every grid cell gets its cap, and the unit total is not reached.
#' As long as the grid cells with a prior can hold the total, this is the solution of
#' \code{\link{toolAllocateByCrossEntropy}} with one hard cap per grid cell and column, computed directly. The
#' spill-over to grid cells without a prior has no counterpart in \code{\link{toolAllocateByCrossEntropy}}, where a
#' zero prior stays zero; use \code{\link{toolAllocateByCrossEntropy}} for caps shared by several flows, soft or
#' signed caps.
#'
#' @param totalLr Numeric matrix of unit totals (units x columns) with the unit names as row names.
#' @param priorHr Numeric matrix of prior weights (grid cells x columns) with grid cell names as row
#' names and the column names of \code{totalLr}.
#' @param capHr Numeric matrix of caps with the dimensions and names of \code{priorHr}. Grid cells with a
#' prior or a cap of zero or below (for example -1e-17 from rounding) take nothing.
#' @inheritParams toolAllocateByCrossEntropy
#' @return Numeric matrix of the grid cell values (grid cells x columns) with the names of \code{priorHr}.
#' @author Patrick v. Jeetze
#' @seealso \code{\link{toolAllocateByCrossEntropy}}
#' @examples
#' total <- matrix(1.5, 1, 1, dimnames = list("A", "protected"))
#' prior <- matrix(c(0.2, 0.1, 0), 3, 1, dimnames = list(c("c1", "c2", "c3"), "protected"))
#' cap <- matrix(c(1, 1, 1), 3, 1, dimnames = dimnames(prior))
#' mapping <- data.frame(cell = c("c1", "c2", "c3"), cluster = "A")
#' # c1 and c2 are filled first (c1 reaches its cap), the rest goes to c3
#' toolAllocateByFilling(total, prior, cap, mapping)
#' @export
toolAllocateByFilling <- function(totalLr, priorHr, capHr, mapping) {
  hrToLr <- as.character(mapping[["cluster"]][match(rownames(priorHr), mapping[["cell"]])])
  stopifnot("priorHr needs the grid cell names of mapping as row names" = !anyNA(hrToLr),
            "totalLr needs all units of the grid cells as row names" = all(hrToLr %in% rownames(totalLr)),
            "priorHr and capHr need the column names of totalLr" =
              all(colnames(totalLr) %in% colnames(priorHr)) & identical(dimnames(capHr), dimnames(priorHr)))

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

  out <- priorHr[, colnames(totalLr), drop = FALSE] * 0
  cellsOf <- split(seq_along(hrToLr), hrToLr)
  for (column in colnames(totalLr)) {
    for (unit in names(cellsOf)) {
      cells <- cellsOf[[unit]]
      total <- totalLr[unit, column]
      x <- .fillCells(priorHr[cells, column], capHr[cells, column], total)
      rest <- total - sum(x)
      if (rest > 0) {
        room <- capHr[cells, column] - x
        x <- x + .fillCells(room, room, rest)
      }
      out[cells, column] <- x
    }
  }
  return(out)
}
