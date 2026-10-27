#' Allocate by cross-entropy
#'
#' Distributes totals of coarse units (Lr, for example the clusters of a model) to their grid cells
#' (Hr), for any quantities given as columns ("flows"), for example land-use transitions or the area
#' of peatland categories. The grid cell values
#' \itemize{
#'   \item add up exactly to the unit totals,
#'   \item stay below upper bounds (caps) on sums of groups of flows in each grid cell, and
#'   \item are as close as possible to a prior distribution (minimum cross-entropy).
#' }
#' The solution is found by iterative proportional fitting extended to upper bounds (see details).
#' For one cap per grid cell and flow, \code{\link{toolAllocateByFilling}} gives the same solution in
#' closed form.
#'
#' @details
#' Minimum cross-entropy means that the function minimises
#' \deqn{\sum_{c,n} X_{cn} \ln(X_{cn} / X^0_{cn}) - X_{cn}}{sum over grid cells c and flows n of
#' X[c,n] * ln(X[c,n] / prior[c,n]) - X[c,n]}
#' subject to the unit totals and the caps.
#'
#' The solution is the prior multiplied by one factor per unit and total and by one factor per
#' grid cell and cap group. The function finds these factors by repeating two scaling steps:
#' \enumerate{
#'   \item For each cap block, groups of flows above their cap are scaled to the cap: flows with
#'   sign 1 by a factor below 1, flows with sign -1 by its inverse. For groups below their cap,
#'   earlier scaling is reversed up to the cap. The accumulated factor of each grid cell and group
#'   (\code{capFactor}) therefore stays between \code{minFactor} and 1.
#'   \item The flows are scaled to the unit totals.
#' }
#' The repetition stops when step 1 no longer changes the unit sums by more than \code{tol}.
#' Because the last step is always the scaling to the unit totals, these hold exactly (unless
#' all grid cells of a unit have a zero cap for a flow); the caps hold up to \code{tol}. The method is
#' iterative proportional fitting extended to upper bounds (Bregman's method). The optimum is
#' unique, so the result does not depend on the order of the cap blocks.
#'
#' A cap block with \code{minFactor} above 0 is soft: where its caps cannot hold, the factor stops
#' at \code{minFactor}, and the cap is exceeded where the objective increases least. Area beyond the
#' cap enters the objective like a flow with a prior smaller by \code{minFactor} that is free of the
#' cap (generalised cross-entropy with error terms). Units where hard caps (\code{minFactor} 0) cannot hold do
#' not converge and show a residual of at least \code{tol}.
#'
#' Only the relative prior within a unit and total matters: multiplying the prior of all grid cells
#' of a unit for one total by a constant does not change the result.
#'
#' @param totalLr magclass object of unit totals (units x one year x flows).
#' @param priorHr magclass object of non-negative prior weights (grid cells x the year of \code{totalLr} x
#' the flows of \code{totalLr}).
#' @param mapping Data frame of the grid cells and their coarse units, as the clustermap of a MAgPIE run,
#' with the columns \code{cell} (names of the grid cells) and \code{cluster} (names of the coarse units). The
#' grid cells are matched by name.
#' @param caps List of cap blocks. Each block is a list with \code{group}, a character vector with
#' one entry per flow naming its group (\code{NA} for flows without cap in this block), \code{cap},
#' a magclass object (the grid cells of \code{priorHr} x the year x groups), and optionally
#' \code{sign}, a vector with one entry per flow (1 or -1, default 1), and \code{minFactor}, the
#' lowest value of the factor of the block (default 0, a hard cap). The sum of a group is the sum
#' of its flows times their sign, so a cap can limit a net change (inflows minus outflows); such a
#' cap can be negative. Groups of one block are disjoint; overlapping constraints (for example
#' source and destination caps of transitions) go into separate blocks.
#' @param fallbackHr magclass object of non-negative weights with the dimensions of \code{priorHr}, used
#' where the prior of a unit and total sums to zero although the total is positive. \code{NULL}
#' gives an error in that case.
#' @param tol Tolerance for the convergence of the caps (absolute, in the measure of \code{totalLr}).
#' @param maxIter Maximum number of iterations.
#' @return List with
#' \item{flowsHr}{magclass object of grid cell flows with the dimensions of \code{priorHr}}
#' \item{capFactor}{list with one magclass object per cap block (dimensions of its \code{cap}): the
#' accumulated factor of each cap group, between \code{minFactor} and 1; below 1 where the cap binds}
#' \item{residual}{magclass object with, per unit, the maximum deviation from the unit totals before the last
#' scaling to them; below \code{tol} where the solution converged}
#' \item{iterations}{number of iterations}
#' @author Patrick v. Jeetze
#' @examples
#' total <- new.magpie("A", NULL, c("crop.past", "crop.forestry"), fill = c(2, 1))
#' prior <- new.magpie(c("c1", "c2", "c3"), NULL, getNames(total), fill = 1)
#' cap <- new.magpie(c("c1", "c2", "c3"), NULL, "crop", fill = c(0.5, 2, 2))
#' source <- list(group = c("crop", "crop"), cap = cap)
#' mapping <- data.frame(cell = c("c1", "c2", "c3"), cluster = "A")
#' toolAllocateByCrossEntropy(total, prior, mapping, caps = list(source))[["flowsHr"]]
#' @export
toolAllocateByCrossEntropy <- function(totalLr, priorHr, mapping, caps = list(), fallbackHr = NULL,
                                       tol = 1e-8, maxIter = 5000) {
  # The iteration runs on matrices (grid cells x flows): written in magclass, it takes about 35 times as long,
  # mostly for selecting grid cells and flows by name in every step.
  .asMatrix <- function(x) {
    return(matrix(x, nrow = dim(x)[1], dimnames = list(getCells(x), getNames(x))))
  }

  # Factor a per grid cell that scales a group above its cap to the cap: the flows with sign 1 (sum plus)
  # by a, the flows with sign -1 (sum minus) by 1 / a, with plus * a - minus / a = cap, for grid cells with
  # flows of the group (plus + minus > 0). For groups below their cap, earlier scaling (capFactor below 1)
  # is reversed up to the cap. The accumulated factor stays at least minFactor.
  .capScale <- function(plus, minus, cap, capFactor, minFactor) {
    # root of plus * a^2 - cap * a - minus = 0, in the form that is exact for cap < 0
    root <- sqrt(cap^2 + 4 * plus * minus)
    factor <- rep(Inf, length(cap))
    capBelowZero <- cap < 0
    factor[capBelowZero] <- 2 * minus[capBelowZero] / (root[capBelowZero] - cap[capBelowZero])
    withPlus <- !capBelowZero & plus > 0
    factor[withPlus] <- (cap[withPlus] + root[withPlus]) / (2 * plus[withPlus])
    # capFactor can reach 0 for hard caps (minFactor 0) where the caps cannot hold
    lower <- if (minFactor > 0) minFactor / capFactor else 0
    factor <- pmax(pmin(factor, 1 / capFactor), lower)
    # groups of hard caps with a zero capFactor whose cap does not bind stay as they are
    factor[is.infinite(factor)] <- 1
    return(factor)
  }

  # Scales the flows of each unit to the unit totals. Totals without any grid cell flow stay zero.
  .scaleToTotalsLr <- function(flows, total, lrIndex) {
    sumLr <- rowsum(flows, lrIndex)
    scale <- total / sumLr
    scale[sumLr == 0] <- 0
    return(flows * scale[lrIndex, , drop = FALSE])
  }

  # Stops with a message if the inputs of toolAllocateByCrossEntropy() do not fit together (the conditions are
  # checked in order, so each can rely on the ones before).
  .checkFlowInput <- function(totalLr, priorHr, hrToLr, caps, fallbackHr) {
    stopifnot(
      "totalLr must be a magclass object with one year" = is.magpie(totalLr) & nyears(totalLr) == 1,
      "totalLr contains NA" = !anyNA(totalLr),
      "totalLr contains negative values" = all(totalLr >= 0),
      "priorHr must be a magclass object with the year of totalLr" =
        is.magpie(priorHr) & identical(getYears(priorHr), getYears(totalLr)),
      "priorHr contains NA" = !anyNA(priorHr),
      "priorHr contains negative values" = all(priorHr >= 0),
      "priorHr must have the flows of totalLr in the same order" = identical(getNames(priorHr), getNames(totalLr)),
      "priorHr needs the grid cell names of mapping" = !anyNA(hrToLr),
      "mapping contains units that are not in totalLr" = all(hrToLr %in% getCells(totalLr)),
      "every unit of totalLr needs at least one grid cell" = all(getCells(totalLr) %in% hrToLr),
      "fallbackHr must be a magclass object without NA" =
        is.null(fallbackHr) | (is.magpie(fallbackHr) & !anyNA(fallbackHr)),
      "fallbackHr contains negative values" = all(fallbackHr >= 0),
      "fallbackHr must have the dimensions and names of priorHr" =
        is.null(fallbackHr) | identical(dimnames(fallbackHr), dimnames(priorHr))
    )
    for (block in caps) {
      sign <- block[["sign"]]
      minFactor <- block[["minFactor"]]
      stopifnot(
        "cap block: group needs one entry per flow" = length(block[["group"]]) == ndata(priorHr),
        "cap block: sign needs one entry per flow, 1 or -1" =
          is.null(sign) | isTRUE(length(sign) == ndata(priorHr) & all(sign %in% c(-1, 1))),
        "cap block: minFactor must be one value in [0, 1)" =
          is.null(minFactor) | isTRUE(length(minFactor) == 1 & all(minFactor >= 0 & minFactor < 1)),
        # caps on net changes can be negative
        "cap block: cap must be a magclass object without NA" = is.magpie(block[["cap"]]) & !anyNA(block[["cap"]]),
        "cap block: cap needs the grid cells of priorHr in its order" =
          identical(getCells(block[["cap"]]), getCells(priorHr)),
        "cap block: every group needs an item in cap" =
          all(block[["group"]][!is.na(block[["group"]])] %in% getNames(block[["cap"]]))
      )
    }
  }


  # coarse unit of each grid cell
  hrToLr <- as.character(mapping[["cluster"]][match(getCells(priorHr), mapping[["cell"]])])
  .checkFlowInput(totalLr, priorHr, hrToLr, caps, fallbackHr)
  total <- .asMatrix(totalLr)
  capMatrix <- lapply(caps, function(block) .asMatrix(block[["cap"]]))

  # row of total for each grid cell; for each cap block its groups (cap column, flows with sign 1 and -1)
  # and its lowest factor
  lrIndex <- match(hrToLr, rownames(total))
  capGroups <- lapply(caps, function(block) {
    column <- match(block[["group"]], getNames(block[["cap"]]))
    sign <- if (is.null(block[["sign"]])) rep(1, ndata(priorHr)) else block[["sign"]]
    return(lapply(unique(column[!is.na(column)]), function(g) {
      list(column = g, plus = which(column == g & sign > 0), minus = which(column == g & sign < 0))
    }))
  })
  minFactor <- lapply(caps, function(block) if (is.null(block[["minFactor"]])) 0 else block[["minFactor"]])

  # start values: the prior, and the fallback for each unit and flow whose prior is zero in all
  # grid cells although the total is positive
  flows <- .asMatrix(priorHr)
  noPrior <- (rowsum(flows, lrIndex) == 0 & total > 0)[lrIndex, , drop = FALSE]
  stopifnot("priorHr is zero in all grid cells of a unit for a flow with positive total and no fallbackHr is given" =
              !any(noPrior) | !is.null(fallbackHr))
  if (any(noPrior)) {
    flows[noPrior] <- .asMatrix(fallbackHr)[noPrior]
  }

  # ---- repeat: scale to the caps, then to the unit totals ------------------------------------

  capFactor <- lapply(capMatrix, function(cap) matrix(1, nrow(cap), ncol(cap), dimnames = dimnames(cap)))
  flows <- .scaleToTotalsLr(flows, total, lrIndex)

  residual <- rep(0, nrow(total))
  # units are independent problems; a unit stops when it has converged
  active <- rep(length(caps) > 0, nrow(total))
  iter <- 0

  while (iter < maxIter && any(active)) {
    iter <- iter + 1
    cells <- which(active[lrIndex])
    activeIndex <- match(lrIndex[cells], which(active))
    activeTotal <- total[active, , drop = FALSE]
    activeFlows <- flows[cells, , drop = FALSE]

    # the groups of a block share no flows, so they are scaled one after another, each only in the grid
    # cells with flows of the group and a factor other than 1 (most caps bind in few grid cells)
    for (b in seq_along(caps)) {
      for (group in capGroups[[b]]) {
        g <- group[["column"]]
        plus <- rowSums(activeFlows[, group[["plus"]], drop = FALSE])
        minus <- rowSums(activeFlows[, group[["minus"]], drop = FALSE])
        withFlows <- which(plus + minus > 0)
        factor <- .capScale(plus[withFlows], minus[withFlows], capMatrix[[b]][cells[withFlows], g],
                            capFactor[[b]][cells[withFlows], g], minFactor[[b]])
        rows <- withFlows[factor != 1]
        factor <- factor[factor != 1]
        if (length(rows) > 0) {
          activeFlows[rows, group[["plus"]]] <- activeFlows[rows, group[["plus"]], drop = FALSE] * factor
          # a zero factor occurs only where the flows with sign -1 are zero
          inverse <- 1 / factor
          inverse[factor == 0] <- 1
          activeFlows[rows, group[["minus"]]] <- activeFlows[rows, group[["minus"]], drop = FALSE] * inverse
          capFactor[[b]][cells[rows], g] <- capFactor[[b]][cells[rows], g] * factor
        }
      }
    }

    # deviation of the unit sums from the totals after scaling to the caps
    residual[active] <- apply(abs(rowsum(activeFlows, activeIndex) - activeTotal), 1, max)
    flows[cells, ] <- .scaleToTotalsLr(activeFlows, activeTotal, activeIndex)
    active[residual < tol] <- FALSE
  }

  # unit sums stay below the totals where all grid cells have a zero cap for a flow
  residual <- pmax(residual, apply(abs(rowsum(flows, lrIndex) - total), 1, max))

  if (max(residual) >= tol) {
    warning("toolAllocateByCrossEntropy did not converge in ", sum(residual >= tol), " of ", nrow(total),
            " units (maximum residual ", signif(max(residual), 3), ")")
  }

  # results with the dimensions and names of the inputs
  flowsHr <- priorHr
  flowsHr[] <- flows
  for (b in seq_along(caps)) {
    factor <- caps[[b]][["cap"]]
    factor[] <- capFactor[[b]]
    capFactor[[b]] <- factor
  }
  residualLr <- new.magpie(getCells(totalLr), getYears(totalLr), fill = residual)
  return(list(flowsHr = flowsHr, capFactor = capFactor, residual = residualLr, iterations = iter))
}
