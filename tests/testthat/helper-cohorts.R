# Values by cohort (stocks or flows) of the grid cells c1, c2, ... (cohortHr) and of the units A, B (cohortLr) as
# magclass objects without year, given grid cell by grid cell (unit by unit).
cohorts <- c("ac0", "ac5", "ac10")
cohortHr <- function(...) {
  values <- matrix(c(...), ncol = 3, byrow = TRUE)
  return(new.magpie(paste0("c", seq_len(nrow(values))), NULL, cohorts, fill = values))
}
cohortLr <- function(...) {
  values <- matrix(c(...), ncol = 3, byrow = TRUE)
  return(new.magpie(c("A", "B")[seq_len(nrow(values))], NULL, cohorts, fill = values))
}
