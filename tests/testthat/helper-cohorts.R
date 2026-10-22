# Cohort stocks of the grid cells c1, c2, ... and unit flows of the units A, B as magclass objects without
# year, values given grid cell by grid cell (unit by unit).
cohorts <- c("ac0", "ac5", "ac10")
cohortStock <- function(...) {
  values <- matrix(c(...), ncol = 3, byrow = TRUE)
  return(new.magpie(paste0("c", seq_len(nrow(values))), NULL, cohorts, fill = values))
}
cohortFlow <- function(...) {
  values <- matrix(c(...), ncol = 3, byrow = TRUE)
  return(new.magpie(c("A", "B")[seq_len(nrow(values))], NULL, cohorts, fill = values))
}
