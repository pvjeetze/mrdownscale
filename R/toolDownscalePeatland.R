#' Downscale MAgPIE peatland
#'
#' Disaggregates the peatland categories of a MAgPIE run (module 58) from clusters to grid cells on a
#' given grid cell land.
#'
#' @details
#' In each year, the categories follow the cluster values: each grid cell keeps its peat area, categories
#' that shrink in the cluster shrink in each grid cell, and the freed area goes to the categories that
#' grow (as in \code{\link[luscale]{interpolate2}}). Where \code{rewettingHr} gives a land reduction on
#' drained peat in a year, the drained peat of each grid cell shrinks by the reduction placed there
#' instead. The result is fitted to the cluster values with
#' \code{\link{toolAllocateByCrossEntropy}}, with the peat area of
#' each grid cell as cap and the drained peat under cropland, pasture and forestry within the grid cell
#' area of that land type (soft cap).
#'
#' Each year depends only on the peat before it, so \code{\link{calcLandTransitionsHighRes}} calls this
#' function year by year together with the land, and the same call on its outputs gives the same
#' peatland.
#'
#' @param landHr magclass object of grid cell land (Mha) by land type, with the years of
#' \code{peatlandLr}. Sub-types (for example "forestry.aff" in \code{landHr} of
#' \code{\link{calcLandTransitionsHighRes}}) are summed to their land type.
#' @param peatlandLr magclass object with the peatland area of each cluster, year and category
#' (as from \code{magpie4::PeatlandArea(gdx, level = "cell", sum = FALSE)}, which is cluster level).
#' @param peatlandStartHr magclass object with the peatland area of each grid cell by category before
#' the first year; missing categories are set to 0.
#' @inheritParams toolAllocateByCrossEntropy
#' @param rewettingHr Optional magclass object with the land reduction on drained peat of each grid cell
#' and year by land type (crop, past, forestry). From \code{transitionsHr} of
#' \code{\link{calcLandTransitionsHighRes}}: \code{dimSums(transitionsHr[, , "peatRewetting"], dim =
#' c("stage", "fromType", "to", "toType"))}.
#' @return List with
#' \item{peatlandHr}{magclass object of grid cell peatland by category (Mha)}
#' \item{excess}{data frame of drained peat above the grid cell land, by year and category (Mha)}
#' @author Patrick v. Jeetze
#' @seealso \code{\link{calcLandTransitionsHighRes}}
#' @examples
#' peatlandLr <- new.magpie("A.1", "y2000", c("crop", "rewetted"), fill = c(0.5, 0.5))
#' peatlandStartHr <- new.magpie(c("c1", "c2"), "y1995", "crop", fill = 0.5)
#' landHr <- new.magpie(c("c1", "c2"), "y2000", "crop", fill = 1)
#' rewettingHr <- new.magpie(c("c1", "c2"), "y2000", "crop", fill = c(0, 0.5))
#' # the cropland on peat of grid cell c2 is reduced, so its drained peat is rewetted
#' mapping <- data.frame(cell = c("c1", "c2"), cluster = "A.1")
#' toolDownscalePeatland(landHr, peatlandLr, peatlandStartHr, mapping, rewettingHr)
#' @export
toolDownscalePeatland <- function(landHr, peatlandLr, peatlandStartHr, mapping, rewettingHr = NULL) {
  cells <- getCells(landHr)
  cluster <- as.character(mapping[["cluster"]][match(cells, mapping[["cell"]])])
  stopifnot("mapping lacks grid cells of landHr" = !anyNA(cluster))
  years <- getYears(peatlandLr)
  categories <- getNames(peatlandLr)
  drained <- intersect(categories, c("crop", "past", "forestry"))
  landType <- getItems(landHr, dim = 3.1, full = TRUE)

  # Denominator of a share: values that are not positive become 1.
  .denominator <- function(x) {
    x[x <= 0] <- 1
    return(x)
  }

  # land of the drained categories (crop, past, forestry), sub-types summed
  land <- new.magpie(cells, years, drained, fill = 0)
  for (k in drained) {
    land[, , k] <- dimSums(landHr[, years, landType == k], dim = 3)
  }
  # peat of the grid cells before the first year, and the land reduction on drained peat; magclass inputs
  # matched by grid cell, missing categories 0
  current <- new.magpie(cells, NULL, categories, fill = 0)
  for (k in intersect(categories, getNames(peatlandStartHr))) {
    current[, , k] <- setYears(peatlandStartHr[cells, 1, k], NULL)
  }
  loss <- land * 0
  if (!is.null(rewettingHr)) {
    for (k in intersect(drained, getNames(rewettingHr))) {
      loss[, , k] <- rewettingHr[cells, years, k]
    }
  }

  out <- list()
  excess <- list()
  for (year in years) {
    # grid cells with peat, and all grid cells of clusters without such grid cells
    hasPeat <- as.vector(dimSums(current, dim = 3)) > 0
    active <- hasPeat | !(cluster %in% cluster[hasPeat])
    activeMapping <- data.frame(cell = cells[active], cluster = cluster[active])
    .onCells <- function(xLr) {
      return(toolAggregate(xLr, activeMapping, from = "cluster", to = "cell"))
    }
    peat <- current[active, , ]
    total <- setYears(peatlandLr[, year, categories], NULL)

    # one step of interpolate2 as prior: in each grid cell, categories shrink by the share by which they
    # shrink in the cluster (drained peat by the land reduction on it, where given), and the freed
    # area goes to the growing categories in proportion to their cluster growth
    previousTotal <- toolAggregate(peat, activeMapping, from = "cell", to = "cluster")[getCells(total), , ]
    change <- total - previousTotal
    shrinkShare <- pmax(-change, 0) / .denominator(previousTotal)
    growShare <- pmax(change, 0) / .denominator(dimSums(pmax(-change, 0), dim = 3))
    shrink <- peat * .onCells(shrinkShare)
    lossYear <- setYears(loss[active, year, ], NULL)
    if (sum(lossYear) > 0) {
      shrink[, , drained] <- pmin(lossYear, peat[, , drained])
    }
    prior <- peat - shrink + dimSums(shrink, dim = 3) * .onCells(growShare)

    # caps: the peat area of each grid cell, and drained peat within the land of its type (soft)
    area <- dimSums(peat, dim = 3)
    caps <- list(list(group = rep("peat", length(categories)), cap = setNames(area, "peat")))
    if (length(drained) > 0) {
      caps[[2]] <- list(group = ifelse(categories %in% drained, categories, NA),
                        cap = setYears(land[active, year, ], NULL), minFactor = 1e-6)
    }
    result <- suppressWarnings(toolAllocateByCrossEntropy(total, prior, mapping, caps, fallbackHr = prior * 0 + area,
                                                          tol = 1e-6, maxIter = 100))

    current[active, , ] <- result[["flowsHr"]]
    out[[year]] <- setYears(current, year)
    above <- pmax(current[, , drained] - setYears(land[, year, ], NULL), 0)
    excess[[year]] <- data.frame(year = rep(year, length(drained)), category = drained,
                                 area = as.vector(dimSums(above, dim = 1)))
  }

  peatlandHr <- mbind(out)
  getSets(peatlandHr, fulldim = FALSE) <- c("x.y.iso", "year", "data")
  return(list(peatlandHr = peatlandHr, excess = do.call(rbind, excess)))
}
