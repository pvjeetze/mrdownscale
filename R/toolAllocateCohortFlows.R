#' Allocate cohort flows to grid cells
#'
#' Distributes flows of cohorts (for example age classes of a land type) of coarse units (Lr, for
#' example the clusters of a model) to their grid cells (Hr), taking them from the grid cell stock of the
#' cohort.
#'
#' @details
#' Without \code{totalHr}, the unit flow of each cohort is distributed in proportion to the
#' grid cell stock of the cohort.
#'
#' With \code{totalHr} (for example the outflow of a land type from the grid cell transitions), the
#' flows of each grid cell add up to its total. The unit flows of the cohorts are first scaled to the
#' sum of the grid cell totals of the unit (where they are all zero, the cohort composition of the
#' stock in the grid cells with a total is used). Starting from the grid cell stock, two steps are repeated (iterative
#' proportional fitting):
#' \enumerate{
#'   \item the flows of each cohort are scaled to its unit flow;
#'   \item the flows of each grid cell are scaled to its total, with each flow at most the stock of its
#'   cohort (flows at the stock stay there, the others take the rest in proportion).
#' }
#' The repetition stops when the unit flows of the cohorts hold up to \code{tol} or no longer
#' improve. Because the grid cell step is the last, the grid cell totals and stocks hold exactly; where the
#' stocks do not allow both, the unit flows of the cohorts deviate. Where the cohorts with a
#' unit flow lack stock in a grid cell, the rest of its total comes from its other cohorts in
#' proportion to their stock.
#'
#' The flow of a grid cell and cohort never exceeds its stock. A unit flow above the unit stock of
#' a cohort, and a grid cell total above the grid cell stock, are not allocated.
#'
#' @param stockHr magclass object of the stock of each grid cell (grid cells x cohorts), without year.
#' @param flowLr magclass object of the flow of each unit (units x the cohorts of \code{stockHr}), without year.
#' @inheritParams toolAllocateByCrossEntropy
#' @param totalHr Optional magclass object with the total flow of each grid cell (the grid cells of
#' \code{stockHr} in its order, one item).
#' @param tol Tolerance for the unit flows of the cohorts (absolute, in the measure of
#' \code{flowLr}).
#' @param maxIter Maximum number of iterations.
#' @return magclass object of the grid cell flows with the dimensions of \code{stockHr}.
#' @author Patrick v. Jeetze
#' @seealso \code{\link{toolShiftCohorts}}, \code{\link{calcAgeClassesHighRes}}
#' @examples
#' stockHr <- new.magpie(c("c1", "c2"), NULL, c("ac0", "ac5"), fill = c(1, 1, 0, 1))
#' flowLr <- new.magpie("A", NULL, c("ac0", "ac5"), fill = c(0.5, 0.5))
#' totalHr <- new.magpie(c("c1", "c2"), NULL, fill = c(0.5, 0.5))
#' # grid cell 1 has no stock of ac5, so its total of 0.5 must come from ac0
#' toolAllocateCohortFlows(stockHr, flowLr, data.frame(cell = c("c1", "c2"), cluster = "A"), totalHr)
#' @export
toolAllocateCohortFlows <- function(stockHr, flowLr, mapping, totalHr = NULL, tol = 1e-8, maxIter = 100) {
  # Scales each row of x to its total, with each entry at most its cap (water filling): entries that
  # the common factor of the row would lift above their cap stay at the cap, and the other entries
  # take the rest of the total in proportion. A total above the sum of the caps is not reached.
  .fillRows <- function(x, cap, total) {
    xSum <- rowSums(x)
    out <- (x / replace(xSum, xSum <= 0, 1)) * total
    atCap <- out > cap
    rows <- which(rowSums(atCap) > 0)
    # each round puts at least one more entry of the remaining rows at its cap
    while (length(rows) > 0) {
      capped <- atCap[rows, , drop = FALSE]
      free <- x[rows, , drop = FALSE] * !capped
      rest <- pmax(total[rows] - rowSums(cap[rows, , drop = FALSE] * capped), 0)
      freeSum <- rowSums(free)
      filled <- (free / replace(freeSum, freeSum <= 0, 1)) * rest + cap[rows, , drop = FALSE] * capped
      above <- filled > cap[rows, , drop = FALSE] & !capped
      out[rows, ] <- filled
      atCap[rows, ] <- capped | above
      rows <- rows[rowSums(above) > 0]
    }
    return(pmin(out, cap))
  }

  # The fit runs on matrices (grid cells x cohorts), as the iteration of toolAllocateByCrossEntropy, which
  # takes about 35 times as long in magclass.
  .asMatrix <- function(x) {
    return(matrix(x, nrow = dim(x)[1], dimnames = list(getCells(x), getNames(x))))
  }
  # result with the dimensions and names of stockHr
  .asStock <- function(x) {
    out <- stockHr
    out[] <- x
    return(out)
  }

  hrToLr <- as.character(mapping[["cluster"]][match(getCells(stockHr), mapping[["cell"]])])
  stopifnot("stockHr and flowLr must be magclass objects" = is.magpie(stockHr) & is.magpie(flowLr),
            "stockHr needs the grid cell names of mapping" = !anyNA(hrToLr),
            "stockHr and flowLr must have the same cohorts" = identical(getNames(stockHr), getNames(flowLr)),
            "totalHr needs the grid cells of stockHr in its order" =
              is.null(totalHr) || identical(getCells(totalHr), getCells(stockHr)))
  stockAll <- .asMatrix(stockHr)
  flowUnits <- .asMatrix(flowLr)

  if (is.null(totalHr)) {
    # in proportion to the grid cell stock of each cohort
    stockOfUnit <- rowsum(stockAll, hrToLr)[rownames(flowUnits), , drop = FALSE]
    share <- pmin(flowUnits, stockOfUnit) / replace(stockOfUnit, stockOfUnit <= 0, 1)
    return(.asStock(stockAll * share[hrToLr, , drop = FALSE]))
  }

  flow <- stockAll * 0
  totalAll <- as.vector(totalHr)
  active <- totalAll > 0
  if (!any(active)) {
    return(.asStock(flow))
  }
  stock <- stockAll[active, , drop = FALSE]
  total <- totalAll[active]
  hrToLr <- hrToLr[active]

  # unit flows of the cohorts scaled to the grid cell totals, within the stock of the active grid cells
  totalSum <- tapply(total, hrToLr, sum)
  unitsLr <- names(totalSum)
  activeStock <- rowsum(stock, hrToLr)[unitsLr, , drop = FALSE]
  target <- flowUnits[unitsLr, , drop = FALSE]
  noFlow <- rowSums(target) == 0
  target[noFlow, ] <- activeStock[noFlow, ]
  targetSum <- rowSums(target)
  target <- pmin(target * (as.vector(totalSum) / replace(targetSum, targetSum <= 0, 1)), activeStock)

  # cohorts with stock in the active grid cells
  cohorts <- colSums(activeStock) > 0
  fitStock <- stock[, cohorts, drop = FALSE]
  target <- target[, cohorts, drop = FALSE]
  x <- .fillRows(fitStock, fitStock, total)
  deviation <- Inf
  for (iter in seq_len(maxIter)) {
    sumLr <- rowsum(x, hrToLr)[unitsLr, , drop = FALSE]
    previous <- deviation
    deviation <- max(0, abs(sumLr - target))
    if (min(deviation, previous - deviation) < tol) {
      break
    }
    sumHr <- sumLr[hrToLr, , drop = FALSE]
    x <- (x / replace(sumHr, sumHr <= 0, 1)) * target[hrToLr, , drop = FALSE]
    x <- .fillRows(x, fitStock, total)
  }

  # rest of the grid cell totals from the other cohorts of the grid cell, in proportion to their stock
  missing <- pmax(total - rowSums(x), 0)
  remaining <- fitStock - x
  remainingSum <- rowSums(remaining)
  x <- x + remaining * pmin(missing / replace(remainingSum, remainingSum <= 0, 1), 1)

  # within the stock also after rounding, so that the stock minus the flow is never negative
  flow[which(active), cohorts] <- pmin(x, fitStock)
  return(.asStock(flow))
}
