#' Shift cohorts
#'
#' Moves the stock of each cohort \code{nShift} cohorts on (for example age classes after a time
#' step); the last cohort keeps its stock and takes the stock shifted beyond it.
#'
#' @param stockHr magclass object of the stock of each grid cell (grid cells x cohorts), cohorts ordered from
#' young to old.
#' @param nShift Number of cohorts to shift (a whole number of at least 0).
#' @return magclass object of the shifted stock; the first \code{nShift} cohorts are empty.
#' @author Patrick v. Jeetze
#' @seealso \code{\link{toolAllocateCohortFlows}}
#' @examples
#' toolShiftCohorts(new.magpie("c1", NULL, c("ac0", "ac5", "acx"), fill = c(1, 2, 3)), 1)
#' @export
toolShiftCohorts <- function(stockHr, nShift) {
  if (nShift == 0) {
    return(stockHr)
  }
  cohorts <- getNames(stockHr)
  n <- length(cohorts)
  shifted <- stockHr * 0
  if (nShift < n) {
    # each cohort takes the stock of the cohort nShift places younger
    moved <- stockHr[, , 1:(n - nShift)]
    getNames(moved) <- cohorts[(nShift + 1):n]
    shifted[, , getNames(moved)] <- moved
  }
  # the last cohort also keeps its stock and takes the stock shifted beyond it
  shifted[, , n] <- shifted[, , n] + dimSums(stockHr[, , max(1, n - nShift + 1):n], dim = 3)
  return(shifted)
}
