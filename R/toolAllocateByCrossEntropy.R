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
#' @param totalLr Numeric matrix of unit totals (units x flows) with row names (units) and
#' column names.
#' @param priorHr Numeric matrix of non-negative prior weights (grid cells x flows) with grid cell names as row
#' names and the column names of \code{totalLr}.
#' @param mapping Data frame of the grid cells and their coarse units, as the clustermap of a MAgPIE run,
#' with the columns \code{cell} (names of the grid cells) and \code{cluster} (names of the coarse units). The
#' grid cells are matched by name (row names of a matrix, cell names of a magclass object).
#' @param caps List of cap blocks. Each block is a list with \code{group}, a character vector with
#' one entry per flow naming its group (\code{NA} for flows without cap in this block), \code{cap},
#' a numeric matrix (grid cells x groups) with the group names as column names, and optionally
#' \code{sign}, a vector with one entry per flow (1 or -1, default 1), and \code{minFactor}, the
#' lowest value of the factor of the block (default 0, a hard cap). The sum of a group is the sum
#' of its flows times their sign, so a cap can limit a net change (inflows minus outflows); such a
#' cap can be negative. Groups of one block are disjoint; overlapping constraints (for example
#' source and destination caps of transitions) go into separate blocks.
#' @param fallbackHr Numeric matrix of non-negative weights with the dimensions of \code{priorHr}, used
#' where the prior of a unit and total sums to zero although the total is positive. \code{NULL}
#' gives an error in that case.
#' @param tol Tolerance for the convergence of the caps (absolute, in the measure of \code{totalLr}).
#' @param maxIter Maximum number of iterations.
#' @return List with
#' \item{flowsHr}{matrix of grid cell flows (grid cells x flows)}
#' \item{capFactor}{list with one matrix (grid cells x groups) per cap block: the accumulated factor of
#' each cap group, between \code{minFactor} and 1; below 1 where the cap binds}
#' \item{residual}{per unit, the maximum deviation from the unit totals before the last
#' scaling to them; below \code{tol} where the solution converged}
#' \item{iterations}{number of iterations}
#' @author Patrick v. Jeetze
#' @examples
#' total <- matrix(c(2, 1), nrow = 1, dimnames = list("A", c("crop.past", "crop.forestry")))
#' prior <- matrix(1, nrow = 3, ncol = 2, dimnames = list(c("c1", "c2", "c3"), colnames(total)))
#' source <- list(group = c("crop", "crop"),
#'                cap = matrix(c(0.5, 2, 2), ncol = 1, dimnames = list(NULL, "crop")))
#' mapping <- data.frame(cell = c("c1", "c2", "c3"), cluster = "A")
#' toolAllocateByCrossEntropy(total, prior, mapping, caps = list(source))[["flowsHr"]]
#' @export
toolAllocateByCrossEntropy <- function(totalLr, priorHr, mapping, caps = list(), fallbackHr = NULL,
                                       tol = 1e-8, maxIter = 5000) {
  # Scales groups above their cap to the cap: flows with sign 1 by a factor a, flows with sign -1 by
  # 1 / a, with plus * a - minus / a = cap. For groups below their cap, earlier scaling (capFactor
  # below 1) is reversed up to the cap. The accumulated factor stays at least minFactor.
  .scaleToCaps <- function(flows, cap, capGroup, capSign, capFactor, minFactor) {
    # sums of the flows with sign 1 (plus) and -1 (minus) in each group (grid cells x groups)
    plus <- matrix(0, nrow = nrow(flows), ncol = ncol(cap))
    minus <- plus
    for (g in unique(capGroup[!is.na(capGroup)])) {
      plus[, g] <- rowSums(flows[, which(capGroup == g & capSign > 0), drop = FALSE])
      minus[, g] <- rowSums(flows[, which(capGroup == g & capSign < 0), drop = FALSE])
    }

    # root of plus * a^2 - cap * a - minus = 0, in the form that is exact for cap < 0
    root <- sqrt(cap^2 + 4 * plus * minus)
    factor <- ifelse(cap < 0, 2 * minus / (root - cap), ifelse(plus > 0, (cap + root) / (2 * plus), Inf))
    # capFactor can reach 0 for hard caps (minFactor 0) where the caps cannot hold
    lower <- if (minFactor > 0) minFactor / capFactor else 0
    factor <- pmax(pmin(factor, 1 / capFactor), lower)
    # groups without flow, and groups of hard caps with a zero capFactor whose cap does not bind,
    # stay as they are
    factor[plus + minus == 0 | is.infinite(factor)] <- 1

    capped <- which(!is.na(capGroup))
    scale <- factor[, capGroup[capped], drop = FALSE]
    # a zero factor occurs only where the flows with sign -1 are zero
    negative <- capSign[capped] < 0
    scale[, negative] <- ifelse(scale[, negative] > 0, 1 / scale[, negative], 1)
    flows[, capped] <- flows[, capped, drop = FALSE] * scale
    return(list(flows = flows, capFactor = capFactor * factor))
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
  .checkFlowInput <- function(total, priorHr, hrToLr, caps, fallbackHr) {
    stopifnot(
      "totalLr must be a numeric matrix" = is.matrix(total) & is.numeric(total),
      "totalLr contains NA" = !anyNA(total),
      "totalLr contains negative values" = all(total >= 0),
      "totalLr needs row names (units) and column names (flows)" =
        !is.null(rownames(total)) & !is.null(colnames(total)),
      "priorHr must be a numeric matrix" = is.matrix(priorHr) & is.numeric(priorHr),
      "priorHr contains NA" = !anyNA(priorHr),
      "priorHr contains negative values" = all(priorHr >= 0),
      "priorHr must have the column names of totalLr in the same order" = identical(colnames(priorHr), colnames(total)),
      "priorHr needs the grid cell names of mapping as row names" = length(hrToLr) == nrow(priorHr) & !anyNA(hrToLr),
      "mapping contains units that are not row names of totalLr" = all(hrToLr %in% rownames(total)),
      "every unit of totalLr needs at least one grid cell" = all(rownames(total) %in% hrToLr),
      "fallbackHr must be a numeric matrix without NA" =
        is.null(fallbackHr) | isTRUE(is.matrix(fallbackHr) & is.numeric(fallbackHr) & !anyNA(fallbackHr)),
      "fallbackHr contains negative values" = all(fallbackHr >= 0),
      "fallbackHr must have the dimensions of priorHr" = is.null(fallbackHr) | identical(dim(fallbackHr), dim(priorHr))
    )
    for (block in caps) {
      sign <- block[["sign"]]
      minFactor <- block[["minFactor"]]
      stopifnot(
        "cap block: group needs one entry per flow" = length(block[["group"]]) == ncol(priorHr),
        "cap block: sign needs one entry per flow, 1 or -1" =
          is.null(sign) | isTRUE(length(sign) == ncol(priorHr) & all(sign %in% c(-1, 1))),
        "cap block: minFactor must be one value in [0, 1)" =
          is.null(minFactor) | isTRUE(length(minFactor) == 1 & all(minFactor >= 0 & minFactor < 1)),
        # caps on net changes can be negative
        "cap block: cap must be a numeric matrix without NA" =
          is.matrix(block[["cap"]]) & is.numeric(block[["cap"]]) & !anyNA(block[["cap"]]),
        "cap block: cap needs one row per grid cell" = nrow(block[["cap"]]) == nrow(priorHr),
        "cap block: every group needs a column in cap" =
          all(block[["group"]][!is.na(block[["group"]])] %in% colnames(block[["cap"]]))
      )
    }
  }


  # coarse unit of each grid cell (row of prior)
  hrToLr <- as.character(mapping[["cluster"]][match(rownames(priorHr), mapping[["cell"]])])
  .checkFlowInput(totalLr, priorHr, hrToLr, caps, fallbackHr)

  # row of totalLr for each grid cell, cap column, sign and lowest factor for each flow
  lrIndex <- match(hrToLr, rownames(totalLr))
  capGroup <- lapply(caps, function(block) match(block[["group"]], colnames(block[["cap"]])))
  capSign <- lapply(caps, function(block) if (is.null(block[["sign"]])) rep(1, ncol(priorHr)) else block[["sign"]])
  minFactor <- lapply(caps, function(block) if (is.null(block[["minFactor"]])) 0 else block[["minFactor"]])

  # start values: the prior, and the fallback for each unit and flow whose prior is zero in all
  # grid cells although the total is positive
  flows <- priorHr
  noPrior <- (rowsum(priorHr, lrIndex) == 0 & totalLr > 0)[lrIndex, , drop = FALSE]
  stopifnot("priorHr is zero in all grid cells of a unit for a flow with positive total and no fallbackHr is given" =
              !any(noPrior) | !is.null(fallbackHr))
  flows[noPrior] <- fallbackHr[noPrior]

  # ---- repeat: scale to the caps, then to the unit totals ------------------------------------

  capFactor <- lapply(caps, function(block) {
    matrix(1, nrow = nrow(block[["cap"]]), ncol = ncol(block[["cap"]]), dimnames = dimnames(block[["cap"]]))
  })
  flows <- .scaleToTotalsLr(flows, totalLr, lrIndex)

  residual <- rep(0, nrow(totalLr))
  # units are independent problems; a unit stops when it has converged
  active <- rep(length(caps) > 0, nrow(totalLr))
  iter <- 0

  while (iter < maxIter && any(active)) {
    iter <- iter + 1
    cells <- active[lrIndex]
    activeIndex <- match(lrIndex[cells], which(active))
    activeTotal <- totalLr[active, , drop = FALSE]
    activeFlows <- flows[cells, , drop = FALSE]

    for (b in seq_along(caps)) {
      scaled <- .scaleToCaps(activeFlows, caps[[b]][["cap"]][cells, , drop = FALSE], capGroup[[b]],
                             capSign[[b]], capFactor[[b]][cells, , drop = FALSE], minFactor[[b]])
      activeFlows <- scaled[["flows"]]
      capFactor[[b]][cells, ] <- scaled[["capFactor"]]
    }

    # deviation of the unit sums from the totals after scaling to the caps
    residual[active] <- apply(abs(rowsum(activeFlows, activeIndex) - activeTotal), 1, max)
    flows[cells, ] <- .scaleToTotalsLr(activeFlows, activeTotal, activeIndex)
    active[residual < tol] <- FALSE
  }

  # unit sums stay below the totals where all grid cells have a zero cap for a flow
  residual <- pmax(residual, apply(abs(rowsum(flows, lrIndex) - totalLr), 1, max))
  names(residual) <- rownames(totalLr)

  if (max(residual) >= tol) {
    warning("toolAllocateByCrossEntropy did not converge in ", sum(residual >= tol), " of ", nrow(totalLr),
            " units (maximum residual ", signif(max(residual), 3), ")")
  }

  return(list(flowsHr = flows, capFactor = capFactor, residual = residual, iterations = iter))
}
