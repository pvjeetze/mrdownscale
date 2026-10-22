# Mapping table of the cells "c1", "c2", ... to the clusters given for each cell.
cellMapping <- function(cluster) {
  return(data.frame(cell = paste0("c", seq_along(cluster)), cluster = cluster))
}
