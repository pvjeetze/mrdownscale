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
#' @param stockHr magclass object of the stock of each grid cell (grid cells x cohorts), without year.
#' @param stockLr magclass object of the stock of each unit (units x the cohorts of \code{stockHr}), without
#' year.
#' @inheritParams toolAllocateByCrossEntropy
#' @param tol Tolerance for the stocks of the units (absolute, in the measure of \code{stockLr}).
#' @param maxIter Maximum number of iterations.
#' @return magclass object of the fitted grid cell stocks with the dimensions of \code{stockHr}.
#' @author Patrick v. Jeetze
#' @seealso \code{\link{toolAllocateCohortFlows}}, \code{\link{calcAgeClassesHighRes}}
#' @examples
#' stockHr <- new.magpie(c("c1", "c2"), NULL, c("ac0", "ac5"), fill = c(1, 0, 0, 1))
#' stockLr <- new.magpie("A", NULL, c("ac0", "ac5"), fill = c(0.5, 1.5))
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
      x[rows, ] <- (x[rows, , drop = FALSE] / replace(rowSum, rowSum <= 0, 1)) * total[rows]
      sumLr <- rowsum(x[rows, , drop = FALSE], lrIndex[rows])
      activeIndex <- which(active)
      deviation <- apply(abs(sumLr - target[activeIndex, , drop = FALSE]), 1, max)
      active[activeIndex[deviation < tol]] <- FALSE
      if (!any(active)) {
        break
      }
      unit <- match(lrIndex[rows], activeIndex)
      sumHr <- sumLr[unit, , drop = FALSE]
      x[rows, ] <- (x[rows, , drop = FALSE] / replace(sumHr, sumHr <= 0, 1)) *
        target[activeIndex[unit], , drop = FALSE]
    }
    return(list(x = x, active = active))
  }

  # The fit runs on matrices (grid cells x cohorts), as the iteration of toolAllocateByCrossEntropy, which
  # takes about 35 times as long in magclass.
  .asMatrix <- function(x) {
    return(matrix(x, nrow = dim(x)[1], dimnames = list(getCells(x), getNames(x))))
  }

  hrToLr <- as.character(mapping[["cluster"]][match(getCells(stockHr), mapping[["cell"]])])
  stopifnot("stockHr and stockLr must be magclass objects" = is.magpie(stockHr) & is.magpie(stockLr),
            "stockHr needs the grid cell names of mapping" = !anyNA(hrToLr),
            "stockHr and stockLr must have the same cohorts" = identical(getNames(stockHr), getNames(stockLr)))
  stock <- .asMatrix(stockHr)
  total <- rowSums(stock)
  totalLr <- tapply(total, hrToLr, sum)
  unitsLr <- names(totalLr)
  lrIndex <- match(hrToLr, unitsLr)

  # cohort composition of the units, scaled to the summed grid cell totals
  target <- .asMatrix(stockLr)[unitsLr, , drop = FALSE]
  targetSum <- rowSums(target)
  share <- target / replace(targetSum, targetSum <= 0, 1)
  share[targetSum == 0, ] <- 1 / ncol(target)
  target <- share * as.vector(totalLr)

  # seed where the fit cannot reach a cohort otherwise: cohorts of a unit that none of its grid cells
  # holds, and grid cells that hold no cohort of their unit
  hasCohort <- target[lrIndex, , drop = FALSE] > 0
  missingLr <- target > 0 & rowsum(stock, lrIndex) == 0
  missingHr <- total > 0 & rowSums(stock * hasCohort) == 0
  seed <- 1e-6 * share[lrIndex, , drop = FALSE] * total
  x <- stock + seed * (missingLr[lrIndex, , drop = FALSE] | missingHr)
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
  out <- stockHr
  out[] <- (fit[["x"]] / replace(fitSum, fitSum <= 0, 1)) * total
  return(out)
}
