#' Fit cohort stocks to unit stocks
#'
#' Fits the stock of each grid cell (Hr) by cohort to the stock of its coarse unit (Lr) by cohort (for
#' example the age classes of a model cluster), keeping the total stock of each grid cell.
#'
#' @details
#' The stocks of the units are scaled to the summed totals of their grid cells (units without
#' stock take equal shares). Where the fit could not reach a cohort otherwise (a cohort of a unit
#' that none of its grid cells holds, or a grid cell that holds no cohort of its unit), a seed of 1e-6 times
#' the cohort composition of the unit times the grid cell total is added; units that do not converge
#' get this seed in all their grid cells and are fitted again. Then two steps
#' are repeated (iterative proportional fitting): the grid cells are scaled to their totals, and the
#' cohorts to the stocks of the units. The repetition stops when the stocks of the units hold up to
#' \code{tol}; the last step scales the grid cells, so the grid cell totals hold exactly.
#'
#' Applied after the cohort flows of a time step (see \code{\link{calcAgeClassesHighRes}}), the
#' area that the fit moves between cohorts, \code{sum(abs(result - stockHr)) / 2}, measures how far
#' the grid cell flows and the stocks of the units disagree.
#'
#' @param stockHr Numeric matrix of the stock of each grid cell (grid cells x cohorts) with grid cell names as row
#' names and cohort names as column names.
#' @param stockLr Numeric matrix of the stock of each unit (units x cohorts) with unit names
#' as row names and the column names of \code{stockHr}.
#' @inheritParams toolAllocateByCrossEntropy
#' @param tol Tolerance for the stocks of the units (absolute, in the measure of \code{stockLr}).
#' @param maxIter Maximum number of iterations.
#' @return Numeric matrix of the fitted grid cell stocks (grid cells x cohorts).
#' @author Patrick v. Jeetze
#' @seealso \code{\link{toolAllocateCohortFlows}}, \code{\link{calcAgeClassesHighRes}}
#' @examples
#' stockHr <- matrix(c(1, 0, 0, 1), 2, 2, dimnames = list(c("c1", "c2"), c("ac0", "ac5")))
#' stockLr <- matrix(c(0.5, 1.5), 1, dimnames = list("A", c("ac0", "ac5")))
#' toolFitCohortStocks(stockHr, stockLr, data.frame(cell = c("c1", "c2"), cluster = "A"))
#' @export
toolFitCohortStocks <- function(stockHr, stockLr, mapping, tol = 1e-8, maxIter = 200) {
  # Iterative proportional fitting of x (grid cells x cohorts) to the grid cell totals (total) and the stocks
  # of the units (target, units x cohorts), for the units marked active; a unit stops when its stocks
  # hold up to tol. Returns x and the units that did not converge (active).
  .fitUnits <- function(x, total, target, lrIndex, active, tol, maxIter) {
    for (iter in seq_len(maxIter)) {
      rows <- which(active[lrIndex])
      rowSum <- rowSums(x[rows, , drop = FALSE])
      x[rows, ] <- (x[rows, , drop = FALSE] / ifelse(rowSum > 0, rowSum, 1)) * total[rows]
      sumLr <- rowsum(x[rows, , drop = FALSE], lrIndex[rows])
      activeIndex <- which(active)
      deviation <- apply(abs(sumLr - target[activeIndex, , drop = FALSE]), 1, max)
      active[activeIndex[deviation < tol]] <- FALSE
      if (!any(active)) {
        break
      }
      unit <- match(lrIndex[rows], activeIndex)
      sumHr <- sumLr[unit, , drop = FALSE]
      x[rows, ] <- (x[rows, , drop = FALSE] / ifelse(sumHr > 0, sumHr, 1)) *
        target[activeIndex[unit], , drop = FALSE]
    }
    return(list(x = x, active = active))
  }

  hrToLr <- as.character(mapping[["cluster"]][match(rownames(stockHr), mapping[["cell"]])])
  stopifnot("stockHr needs the grid cell names of mapping as row names" =
              length(hrToLr) == nrow(stockHr) & !anyNA(hrToLr))
  total <- rowSums(stockHr)
  totalLr <- tapply(total, hrToLr, sum)
  unitsLr <- names(totalLr)
  lrIndex <- match(hrToLr, unitsLr)

  # cohort composition of the units, scaled to the summed grid cell totals
  target <- stockLr[unitsLr, , drop = FALSE]
  targetSum <- rowSums(target)
  share <- target / ifelse(targetSum > 0, targetSum, 1)
  share[targetSum == 0, ] <- 1 / ncol(target)
  target <- share * as.vector(totalLr)

  # seed where the fit cannot reach a cohort otherwise: cohorts of a unit that none of its grid cells
  # holds, and grid cells that hold no cohort of their unit
  hasCohort <- target[lrIndex, , drop = FALSE] > 0
  missingLr <- target > 0 & rowsum(stockHr, lrIndex) == 0
  missingHr <- total > 0 & rowSums(stockHr * hasCohort) == 0
  seed <- 1e-6 * share[lrIndex, , drop = FALSE] * total
  x <- stockHr + seed * (missingLr[lrIndex, , drop = FALSE] | missingHr)
  fit <- .fitUnits(x, total, target, lrIndex, rep(TRUE, length(unitsLr)), tol, maxIter)

  # units that did not converge (grid cells whose cohorts cannot reach the stocks of the unit): again
  # with the seed in all their grid cells
  if (any(fit[["active"]])) {
    rows <- which(fit[["active"]][lrIndex])
    x <- fit[["x"]]
    x[rows, ] <- x[rows, , drop = FALSE] + seed[rows, , drop = FALSE]
    fit <- .fitUnits(x, total, target, lrIndex, fit[["active"]], tol, maxIter)
  }
  fitSum <- rowSums(fit[["x"]])
  return((fit[["x"]] / ifelse(fitSum > 0, fitSum, 1)) * total)
}
