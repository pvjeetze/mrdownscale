# Cohorts and matrices of cohort stocks (grid cells c1, c2, ...) and unit flows (units A, B) by row.
cohorts <- c("ac0", "ac5", "ac10")
cohortStock <- function(...) withCells(matrix(c(...), ncol = 3, byrow = TRUE, dimnames = list(NULL, cohorts)))
cohortFlow <- function(...) {
  return(matrix(c(...), ncol = 3, byrow = TRUE, dimnames = list(c("A", "B")[seq_len(length(c(...)) / 3)], cohorts)))
}
