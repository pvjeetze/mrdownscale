#' Shift cohorts
#'
#' Moves the stock of each cohort \code{nShift} cohorts on (for example age classes after a time
#' step); the last cohort keeps its stock and takes the stock shifted beyond it.
#'
#' @param stockHr Numeric matrix of the stock of each grid cell (grid cells x cohorts), cohorts ordered from
#' young to old.
#' @param nShift Number of cohorts to shift (a whole number of at least 0).
#' @return Numeric matrix of the shifted stock; the first \code{nShift} cohorts are empty.
#' @author Patrick v. Jeetze
#' @seealso \code{\link{toolAllocateCohortFlows}}
#' @examples
#' toolShiftCohorts(matrix(c(1, 2, 3), 1, dimnames = list(NULL, c("ac0", "ac5", "acx"))), 1)
#' @export
toolShiftCohorts <- function(stockHr, nShift) {
  if (nShift == 0) {
    return(stockHr)
  }
  n <- ncol(stockHr)
  shifted <- stockHr * 0
  if (nShift >= n) {
    shifted[, n] <- rowSums(stockHr)
    return(shifted)
  }
  shifted[, (nShift + 1):n] <- stockHr[, 1:(n - nShift), drop = FALSE]
  shifted[, n] <- shifted[, n] + rowSums(stockHr[, (n - nShift + 1):n, drop = FALSE])
  return(shifted)
}
