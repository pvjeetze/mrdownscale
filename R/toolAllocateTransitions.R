#' Allocate land-use transitions
#'
#' Disaggregates the land-use transitions of coarse units (Lr, for example the clusters of a model)
#' to their grid cells (Hr), one time step and stage after the other. The grid cell transitions sum exactly to
#' the unit transitions, so the grid cell states sum to the unit states and the gross change summed over
#' the grid cells of a unit equals the gross change of the unit.
#'
#' @details
#' Each year consists of stages in the order in which the model applies them (for MAgPIE, for
#' example, the shifts before the optimisation and the transitions of the optimisation). Each
#' stage is allocated with \code{\link{toolAllocateByCrossEntropy}}:
#' \itemize{
#'   \item the prior of a transition from type f to type k in grid cell c is
#'   \code{state[c, f] * fromWeight[c, f] * toWeight[c, k]}, so by default the transitions out of a
#'   type are proportional to its area in the grid cell;
#'   \item the area that leaves a type in a grid cell cannot exceed the area of the type in the grid cell;
#'   \item optional \code{constraints} bound, per grid cell, the area of a group of types at the end of
#'   each stage (for example cropland within the suitable area, or protected land). Transitions
#'   into the group count positive, transitions out of it negative, so a grid cell beyond its bound loses
#'   area first. Where the constraints cannot hold together with the unit transitions, they are
#'   exceeded in the grid cells where the objective increases least, constraints with a lower priority
#'   first, and the area that the transitions add beyond them is reported.
#' }
#' Constraints are soft (\code{minFactor} in \code{\link{toolAllocateByCrossEntropy}}): area beyond a constraint
#' enters the objective like a flow with a prior smaller by 1e-6 for the lowest priority, 1e-12 for
#' the next, and so on.
#' The new state of each grid cell is the old state minus the outflows plus the inflows.
#'
#' Before the first year, \code{startHr} is fitted to the unit state before the first year
#' (the unit state at the end of the first year minus the net change of
#' its transitions) by iterative proportional fitting, keeping the area of each grid cell. The largest
#' deviation of the unit sums of the grid cell states from \code{statesLr} is reported.
#'
#' Weights that depend on the state, on the previous year or on input data come from
#' \code{stageRules}, a function with arguments \code{year}, \code{stage}, \code{stateHr} (magclass object,
#' grid cells x the year x types), \code{previousTransitionsHr} (the grid cell transitions of the previous
#' year, a magclass object as \code{transitionsHr} of the return value for that year; \code{NULL} in the first
#' year) and \code{layer} (see below). It returns \code{NULL} or a list with any of \code{fromWeight} and
#' \code{toWeight}, magclass objects (grid cells x types, in the order of the grid cells of \code{startHr});
#' types without an item get weight 1.
#'
#' The list can also contain \code{targetArea}, a magclass object (grid cells x one type): the grid cell
#' distribution of a type at the end of the stage (for example urban land from gridded data), scaled in
#' each unit to the unit area of the type at the end of the stage. The
#' grid cell transitions of this type then follow the target: its gain and loss per grid cell are scaled
#' to the unit inflow and outflow, the gain comes from the types present in the grid cell in proportion
#' to their area, and the
#' other unit transitions are fitted so that the unit states remain those of the model.
#' The transitions that can be added in this fit are given by \code{allowedTransitions}, a character
#' vector of transitions "from.to"; by default these are the transitions with a positive total
#' in the stage. The changes against the unit transitions are reported in
#' \code{transitionChange}.
#'
#' A quantity that follows the land year by year and in turn sets weights of the next year (for
#' example drained peat under managed land) is a \code{layer}: an object with its start value, which
#' \code{updateLayer(year, stateHr, transitionsHr, layer)} updates at the end of each year from
#' the grid cell state and the grid cell transitions of the year (magclass objects as for \code{stageRules}).
#' \code{stageRules} receives the layer of the previous year.
#'
#' @param startHr magclass object of the grid cell state before the first year (grid cells x one year x
#' types); type names must not contain ".".
#' @inheritParams toolAllocateByCrossEntropy
#' @param statesLr magclass object of the unit states at the end of each year (units x the years of
#' \code{transitionsLr} x the types of \code{startHr}), for example the cluster states of a model; for models
#' with states only, the transitions come from \code{\link{toolMinimalTransitions}}.
#' @param transitionsLr magclass object of unit transitions (the units of \code{statesLr} x years x
#' "stage.from.to"), the stages of each year in the order of the model (the order of the items), as from
#' \code{magpie4::landTransitions}. Transitions from a type to itself are ignored.
#' @param stageRules Optional function returning weights and a target area per year and stage, see
#' details.
#' @param constraints Optional named list of constraints. Each is a list with
#' \itemize{
#'   \item \code{types}: the types whose summed grid cell area is bounded;
#'   \item \code{upper} or \code{lower}: the bound, a magclass object (grid cells x years, one item; matched by
#'   grid cell). A year uses the latest year of the bound at or before it, earlier years the first. A lower
#'   bound (for example protected land) acts as an upper bound on all other types: the grid cell area of
#'   \code{startHr} minus the bound;
#'   \item \code{priority}: a number, default 1; constraints with a lower priority are exceeded first.
#' }
#' @param layer Optional start value of a layer, see details.
#' @param updateLayer Optional function that updates the layer at the end of each year, see
#' details.
#' @param tol Tolerance passed to \code{\link{toolAllocateByCrossEntropy}}.
#' @param maxIter Maximum number of iterations passed to \code{\link{toolAllocateByCrossEntropy}}.
#' @return List with
#' \item{startHr}{magclass object of the grid cell state before the first year, fitted to the unit states}
#' \item{stateHr}{magclass object of the grid cell states at the end of each year (grid cells x years x types)}
#' \item{transitionsHr}{magclass object of the grid cell transitions (grid cells x years x "stage.from.to",
#' in the order of the stages; transitions without flow in any grid cell and year are left out)}
#' \item{diagnostics}{data frame by year and stage with the number of iterations, the maximum
#' residual and the number of units where constraints are exceeded}
#' \item{excess}{data frame by year, stage, constraint and unit (column cluster) with the area that the grid cell
#' transitions add beyond the constraints (\code{NULL} where they hold)}
#' \item{transitionChange}{data frame by year, stage and transition with the summed absolute change
#' of the unit transitions due to a target (\code{NULL} without target)}
#' \item{layer}{the layer after the last year (\code{NULL} without layer)}
#' \item{deviationLr}{largest deviation of the unit sums of the grid cell states from \code{statesLr}}
#' @author Patrick v. Jeetze
#' @seealso \code{\link{toolAllocateByCrossEntropy}}
#' @examples
#' startHr <- new.magpie(c("c1", "c2"), "y1995", c("crop", "other"), fill = c(1, 0.5, 0, 0.5))
#' statesLr <- new.magpie("A", "y2000", c("crop", "other"), fill = c(1.1, 0.9))
#' transitionsLr <- new.magpie("A", "y2000", "opt.crop.other", fill = 0.4)
#' mapping <- data.frame(cell = c("c1", "c2"), cluster = "A")
#' toolAllocateTransitions(startHr, mapping, statesLr, transitionsLr)[["stateHr"]]
#' @export
toolAllocateTransitions <- function(startHr, mapping, statesLr, transitionsLr, stageRules = NULL, # nolint: cyclocomp_linter
                                    constraints = NULL, layer = NULL, updateLayer = NULL, tol = 1e-8,
                                    maxIter = 200) {
  # Stops with a message if the inputs of toolAllocateTransitions() do not fit together (the conditions
  # are checked in order, so each can rely on the ones before).
  .checkTransitionInput <- function(startHr, hrToLr, statesLr, transitionsLr, constraints) {
    types <- getNames(startHr)
    stopifnot(
      "startHr must be a magclass object with one year" = is.magpie(startHr) & nyears(startHr) == 1,
      "type names must not contain '.'" = !any(grepl(".", types, fixed = TRUE)),
      "startHr must not contain NA or negative values" = !anyNA(startHr) & all(startHr >= 0, na.rm = TRUE),
      "startHr needs the grid cell names of mapping" = !anyNA(hrToLr),
      "statesLr must be a magclass object with the types of startHr" =
        is.magpie(statesLr) & setequal(getNames(statesLr), types),
      "statesLr must not contain NA values" = !anyNA(statesLr),
      "mapping contains units that are not in statesLr" = all(hrToLr %in% getCells(statesLr)),
      "transitionsLr must be a magclass object with the units and years of statesLr" =
        is.magpie(transitionsLr) & setequal(getCells(transitionsLr), getCells(statesLr)) &
        identical(getYears(transitionsLr), getYears(statesLr)),
      "transitionsLr needs items 'stage.from.to' with types of startHr" =
        all(grepl("^[^.]+[.][^.]+[.][^.]+$", getNames(transitionsLr))) &
        all(getItems(transitionsLr, dim = 3.2, full = TRUE) %in% types) &
        all(getItems(transitionsLr, dim = 3.3, full = TRUE) %in% types),
      "transitionsLr must not contain NA or negative values" =
        !anyNA(transitionsLr) & all(transitionsLr >= 0, na.rm = TRUE),
      "constraints must be a named list" = length(constraints) == 0 | !is.null(names(constraints))
    )
    for (constraint in constraints) {
      bound <- if (is.null(constraint[["upper"]])) constraint[["lower"]] else constraint[["upper"]]
      stopifnot(
        "a constraint contains types not in startHr" = all(constraint[["types"]] %in% types),
        "each constraint needs exactly one of upper and lower" =
          is.null(constraint[["upper"]]) != is.null(constraint[["lower"]]),
        "constraint bounds must be magclass objects with the grid cells of startHr and one item" =
          isTRUE(is.magpie(bound) & ndata(bound) == 1 & all(getCells(startHr) %in% getCells(bound)))
      )
    }
  }

  # The allocation runs on matrices (grid cells or units x types or transitions), converted from and to
  # magclass at the interfaces: in magclass, the kernel iteration takes about 35 times as long (see
  # toolAllocateByCrossEntropy).
  .asMagpie <- function(x, year = NULL) {
    out <- new.magpie(rownames(x), year, colnames(x))
    out[] <- x
    return(out)
  }
  .asMatrix <- function(x) {
    return(matrix(x, nrow = dim(x)[1], dimnames = list(getCells(x), getNames(x))))
  }

  # Unit state after the first nAllocated stages of a year (units x types): the state at the end of the
  # year minus the net change of the later stages; nAllocated = 0 gives the state before the year.
  .unitState <- function(year, nAllocated) {
    stages <- stagesLr[[year]]
    state <- matrix(statesUnits[, , year], nrow = dim(statesUnits)[1], dimnames = dimnames(statesUnits)[1:2])
    for (total in stages[seq_along(stages) > nAllocated]) {
      total <- total[rownames(state), , drop = FALSE]
      from <- sub("[.].*$", "", colnames(total))
      to <- sub("^.*[.]", "", colnames(total))
      for (n in seq_len(ncol(total))) {
        state[, from[n]] <- state[, from[n]] + total[, n]
        state[, to[n]] <- state[, to[n]] - total[, n]
      }
    }

    # rounding in the model output leaves values of about -1e-15
    state[state < 0] <- 0
    return(state)
  }

  # Fits the initial grid cell state to the unit states by iterative proportional fitting: scales each
  # grid cell to its area and each type to its unit total, ending with the unit totals. Stops when
  # the grid cell areas hold or no longer improve (where the unit area of start and the summed grid cell
  # area differ slightly).
  .fitInitial <- function(initial, hrToLr, start, tol, maxIter = 1000) {
    lrIndex <- match(hrToLr, rownames(start))
    start <- start[, colnames(initial), drop = FALSE]
    cellArea <- rowSums(initial)

    # each type scaled to its unit total; types without area in the grid cells of a unit stay zero
    .toUnitTotals <- function(state) {
      sumLr <- rowsum(state, lrIndex)
      scale <- start / sumLr
      scale[sumLr == 0] <- 0
      return(state * scale[lrIndex, , drop = FALSE])
    }

    state <- .toUnitTotals(initial)
    areaDeviation <- max(abs(rowSums(state) - cellArea))
    for (iter in seq_len(maxIter)) {
      state <- .toUnitTotals(state * ifelse(rowSums(state) > 0, cellArea / rowSums(state), 0))
      previousDeviation <- areaDeviation
      areaDeviation <- max(abs(rowSums(state) - cellArea))
      if (min(areaDeviation, previousDeviation - areaDeviation) < tol) {
        break
      }
    }
    return(state)
  }

  # Allocates the unit transitions of one stage between different types with a positive total to
  # grid cells and returns the new grid cell state, the grid cell transitions and the diagnostics: the prior is the
  # area of the source type times the weights of input, the area of the source type is a hard cap,
  # and the constraints of the year (input[["caps"]]) are soft caps.
  .allocateStage <- function(state, mapping, hrToLr, total, input, tol, maxIter) {
    from <- sub("[.].*$", "", colnames(total))
    to <- sub("^.*[.]", "", colnames(total))
    keep <- from != to & colSums(total) > 0

    if (any(keep)) {
      total <- total[, keep, drop = FALSE]
      from <- from[keep]
      to <- to[keep]

      # prior: area of the source type, times the source and destination weights (1 for types
      # without weight)
      fromWeight <- matrix(1, nrow = nrow(state), ncol = ncol(state), dimnames = dimnames(state))
      toWeight <- fromWeight
      fromWeight[, colnames(input[["fromWeight"]])] <- input[["fromWeight"]]
      toWeight[, colnames(input[["toWeight"]])] <- input[["toWeight"]]
      prior <- state[, from, drop = FALSE] * fromWeight[, from, drop = FALSE] * toWeight[, to, drop = FALSE]
      colnames(prior) <- colnames(total)
      # where the weights are zero in all grid cells of a unit, transitions follow the source area
      fallback <- state[, from, drop = FALSE]
      colnames(fallback) <- colnames(total)

      # soft caps of the stage: for each cap of the year whose area the transitions change, the transitions that
      # change it (member), their sign (1 into the types, -1 out of them), the room per grid cell (cap minus the
      # current area of the types, negative where the grid cell is above the cap) and the priority
      softCaps <- list()
      priorities <- numeric(0)
      for (name in names(input[["caps"]])) {
        x <- input[["caps"]][[name]]
        into <- to %in% x[["types"]] & !from %in% x[["types"]]
        outOf <- from %in% x[["types"]] & !to %in% x[["types"]]
        if (any(into | outOf)) {
          room <- as.vector(x[["cap"]]) - rowSums(state[, x[["types"]], drop = FALSE])
          softCaps[[name]] <- list(member = into | outOf, sign = ifelse(outOf, -1, 1), cap = room,
                                   priority = x[["priority"]])
          priorities <- c(priorities, x[["priority"]])
        }
      }
      priorities <- sort(unique(priorities))

      # hard cap: the outflow of a type cannot exceed its area; soft caps with a lower priority get a
      # higher minFactor and are exceeded first
      caps <- list(list(group = from, cap = .asMagpie(state)))
      for (name in names(softCaps)) {
        x <- softCaps[[name]]
        room <- matrix(x[["cap"]], ncol = 1, dimnames = list(rownames(state), name))
        caps[[name]] <- list(group = ifelse(x[["member"]], name, NA), sign = x[["sign"]], cap = .asMagpie(room),
                             minFactor = 1e-6^match(x[["priority"]], priorities))
      }

      # units that do not converge within maxIter are corrected by .fitSourceArea and show in the
      # residual of the diagnostics
      solved <- suppressWarnings(toolAllocateByCrossEntropy(.asMagpie(total), .asMagpie(prior), mapping, caps,
                                                            .asMagpie(fallback), tol, maxIter))
      flows <- matrix(solved[["flowsHr"]], nrow = nrow(state), dimnames = dimnames(prior))
      flows <- .fitSourceArea(flows, state, from, hrToLr, tol)

      result <- list(state = .applyFlows(state, flows), flows = flows, residual = as.vector(solved[["residual"]]),
                     iterations = solved[["iterations"]], excess = .constraintExcess(softCaps, flows, hrToLr, tol))
    } else {
      result <- list(state = state, flows = matrix(0, nrow(state), 0), residual = 0, iterations = 0)
    }
    return(result)
  }

  # Where the outflow of a type exceeds its area in a grid cell (toolAllocateByCrossEntropy stopped before it
  # converged), scales the outflow down to the area and moves the removed area to the grid cells of the
  # unit that have area left, in proportion to that area and at most up to it. The unit totals stay the
  # same unless the unit lacks the room, which is reported.
  .fitSourceArea <- function(flows, state, from, hrToLr, tol) {
    for (type in unique(from)) {
      n <- which(from == type)
      outflow <- rowSums(flows[, n, drop = FALSE])
      over <- outflow > state[, type]
      # outflow above the area of the type
      removed <- flows[, n, drop = FALSE] * (1 - pmin(state[, type] / outflow, 1))
      removed[!over, ] <- 0
      flows[, n] <- flows[, n, drop = FALSE] - removed

      # moved to the grid cells of the unit with area left
      room <- pmax(state[, type] - rowSums(flows[, n, drop = FALSE]), 0)
      room[over] <- 0
      removedCluster <- rowsum(removed, hrToLr)[hrToLr, , drop = FALSE]
      roomCluster <- rowsum(room, hrToLr)[hrToLr, 1]
      removedTotal <- rowSums(removedCluster)
      fill <- ifelse(removedTotal > roomCluster, roomCluster / removedTotal, 1)
      added <- room / replace(roomCluster, roomCluster <= 0, 1) * removedCluster * fill
      flows[, n] <- flows[, n, drop = FALSE] + added

      lost <- sum(removed) - sum(added)
      if (lost > tol) {
        warning("outflow above the area of ", type, " could not be moved: ", signif(lost, 3))
      }
    }
    return(flows)
  }

  # Area that the grid cell transitions add above the soft caps, per constraint and unit (rows above
  # tol).
  .constraintExcess <- function(softCaps, flows, hrToLr, tol) {
    result <- data.frame(constraint = character(0), cluster = character(0), area = numeric(0))
    for (name in names(softCaps)) {
      x <- softCaps[[name]]
      netChange <- as.vector(flows[, x[["member"]], drop = FALSE] %*% x[["sign"]][x[["member"]]])
      excess <- pmax(netChange - pmax(x[["cap"]], 0), 0)
      area <- tapply(excess, hrToLr, sum)
      area <- area[area > tol]
      result <- rbind(result, data.frame(constraint = rep(name, length(area)), cluster = names(area),
                                         area = as.vector(area)))
    }
    rownames(result) <- NULL
    return(result)
  }

  # Grid cell transitions of a type with a target grid cell area, and the other unit transitions balanced
  # so that the unit states stay those of the model. The target distribution of the type (targetHr, a vector
  # by grid cell) is scaled in each unit to the unit area of the type at the end of the stage (unitState), then:
  # 1. the gain and loss per grid cell are the change towards the target, within the area available in the
  #    grid cell; the gain is scaled to keep the model's net change of the unit;
  # 2. the gain comes from the types in the grid cell that can gain or that supply the type in the model, in
  #    proportion to their area; a type that cannot gain gives at most what it gives in the model, the rest of
  #    its share comes from types that can gain;
  # 3. the loss goes to the destination types of the model in proportion to their unit totals, or
  #    to the types in the grid cell that can gain where the model has no outflow;
  # 4. the other transitions are changed by the net extra loss of each type (.rebalance);
  # 5. what they cannot compensate goes back to the target inflow (.swapSources).
  .targetFlows <- function(state, mapping, hrToLr, total, type, targetHr, unitState, allowed, tol) {
    types <- colnames(state)
    unitsLr <- rownames(total)
    lrIndex <- match(hrToLr, unitsLr)

    # target area of each grid cell: its share of the unit target times the unit area of the type
    targetSum <- tapply(targetHr, hrToLr, sum)[unitsLr]
    targetSumHr <- targetSum[hrToLr]
    targetShare <- targetHr / replace(targetSumHr, targetSumHr <= 0, 1)
    target <- as.vector(unitState[hrToLr, type] * targetShare)

    # unit transitions as an array units x from x to, and the names "from.to" of its pairs
    from <- sub("[.].*$", "", colnames(total))
    to <- sub("^.*[.]", "", colnames(total))
    transitions <- array(0, dim = c(nrow(total), length(types), length(types)),
                         dimnames = list(unitsLr, types, types))
    for (n in seq_len(ncol(total))) {
      if (from[n] != to[n]) {
        transitions[, from[n], to[n]] <- transitions[, from[n], to[n]] + total[, n]
      }
    }
    pairs <- paste(rep(types, times = length(types)), rep(types, each = length(types)), sep = ".")
    typePairs <- c(paste(types, type, sep = "."), paste(type, types, sep = "."))
    intoType <- matrix(transitions[, , type], length(unitsLr), dimnames = list(unitsLr, types))
    outOfType <- matrix(transitions[, type, ], length(unitsLr), dimnames = list(unitsLr, types))
    if (is.null(allowed)) {
      allowed <- apply(transitions > 0, c(2, 3), any)
    }
    allowed <- allowed[types, types]
    diag(allowed) <- FALSE
    allowed[type, ] <- FALSE
    allowed[, type] <- FALSE
    # a type can gain if a transition into it is allowed
    canGain <- colSums(allowed) > 0
    canGainCell <- matrix(canGain, nrow(state), length(types), byrow = TRUE)
    canGainUnit <- matrix(canGain, length(unitsLr), length(types), byrow = TRUE)

    # 1. gain and loss per grid cell
    eligible <- (canGainUnit | intoType > 0)
    eligible[, type] <- FALSE
    sourceArea <- state * eligible[lrIndex, , drop = FALSE]
    change <- target - state[, type]
    loss <- pmin(pmax(-change, 0), state[, type])
    # where the model has no outflow of the type, a grid cell loses it only to types in the grid cell that can gain
    noModelOutflow <- (rowSums(outOfType) == 0)[lrIndex]
    gainArea <- state * canGainCell
    gainSum <- rowSums(gainArea)
    loss[noModelOutflow & gainSum == 0] <- 0
    netChange <- rowSums(intoType) - rowSums(outOfType)
    gainTotal <- pmax(netChange + as.vector(rowsum(loss, lrIndex)), 0)
    # unit totals (vector by unit) over the grid cells in proportion to the prior, within the caps
    .byFilling <- function(prior, cap, totalLr) {
      out <- toolAllocateByFilling(.asMagpie(matrix(totalLr, dimnames = list(unitsLr, "x"))),
                                   .asMagpie(matrix(prior, dimnames = list(rownames(state), "x"))),
                                   .asMagpie(matrix(cap, dimnames = list(rownames(state), "x"))), mapping)
      return(as.vector(out))
    }
    gain <- .byFilling(pmax(change, 0), rowSums(sourceArea), gainTotal)

    # 2. sources of the gain
    sourceSum <- rowSums(sourceArea)
    inflow <- gain * (sourceArea / replace(sourceSum, sourceSum <= 0, 1))
    limit <- ifelse(canGainUnit, Inf, intoType)
    supply <- rowsum(inflow, lrIndex)
    supplyScale <- ifelse(supply > 0, pmin(limit / supply, 1), 1)
    if (any(supplyScale < 1)) {
      limited <- inflow * supplyScale[lrIndex, , drop = FALSE]
      excess <- rowSums(inflow) - rowSums(limited)
      otherArea <- (sourceArea - limited) * canGainCell
      moved <- .byFilling(excess, rowSums(otherArea), as.vector(rowsum(excess, lrIndex)))
      otherSum <- rowSums(otherArea)
      inflow <- limited + moved * (otherArea / replace(otherSum, otherSum <= 0, 1))
    }

    # 3. destinations of the loss
    outSum <- rowSums(outOfType)
    destinationShare <- (outOfType / replace(outSum, outSum <= 0, 1))[lrIndex, , drop = FALSE]
    destinationShare[noModelOutflow, ] <- (gainArea / replace(gainSum, gainSum <= 0, 1))[noModelOutflow, ]
    outflow <- loss * destinationShare
    # grid cell transitions into and out of the type ("from.to"), without those that are zero in all grid cells
    flows <- cbind(inflow, outflow)
    colnames(flows) <- typePairs
    flows <- flows[, colSums(flows) > 0, drop = FALSE]

    # 4. other transitions: a type that loses more to the target type, or gains less from it, loses
    #    less elsewhere
    lostMore <- (rowsum(inflow, lrIndex) - intoType) -
      (rowsum(outflow, lrIndex) - outOfType)
    rest <- transitions
    rest[, type, ] <- 0
    rest[, , type] <- 0
    areaAfter <- rowsum(.applyFlows(state, flows), lrIndex)
    balanced <- .rebalance(rest, lostMore, allowed, areaAfter)

    # 5. what the other transitions cannot compensate is taken from the type that has to lose more
    #    instead of the type that has to lose less, in the target inflow of the grid cells
    take <- attr(balanced, "take")
    if (any(take > tol)) {
      inflow <- .swapSources(inflow, state - inflow, take, attr(balanced, "give") * eligible, lrIndex)
      remaining <- attr(inflow, "remaining")
      if (max(remaining) >= tol) {
        warning("transitions could not be balanced in ", sum(remaining >= tol), " units (",
                signif(sum(remaining), 3), " Mha)")
      }
      flows <- cbind(inflow, outflow)
      colnames(flows) <- typePairs
      flows <- flows[, colSums(flows) > 0, drop = FALSE]
    }

    typeChange <- c(colSums(abs(rowsum(inflow, lrIndex) - intoType)),
                    colSums(abs(rowsum(outflow, lrIndex) - outOfType)))
    restChange <- apply(abs(balanced - rest), c(2, 3), sum)
    changeTable <- data.frame(transition = c(typePairs, pairs), change = c(typeChange, as.vector(restChange)))
    changeTable <- changeTable[changeTable[["change"]] > tol, ]

    # the other transitions as a matrix with columns "from.to", for the positive transitions
    rest <- matrix(balanced, nrow = dim(balanced)[1], dimnames = list(unitsLr, pairs))
    return(list(flows = flows, rest = rest[, colSums(rest) > 0, drop = FALSE], change = changeTable))
  }

  # Changes the transitions of each unit so that each type loses lostMore less (lostMore > 0) or
  # -lostMore more (lostMore < 0). Area moves from types that have to lose more (giver) to types that
  # have to lose less (taker), in four steps (.compensate). Returns the transitions with the area that
  # could not be moved (units x types) as attributes "take" and "give".
  .rebalance <- function(transitions, lostMore, allowed, area) {
    take <- pmax(lostMore, 0)
    give <- pmax(-lostMore, 0)
    for (j in which(rowSums(take) > 0)) {
      m <- transitions[j, , ]
      for (step in 1:4) {
        for (taker in which(take[j, ] > 0)) {
          for (giver in which(give[j, ] > 0)) {
            result <- .compensate(step, m, taker, giver, min(take[j, taker], give[j, giver]), allowed, area[j, ])
            m <- result[["m"]]
            take[j, taker] <- take[j, taker] - result[["moved"]]
            give[j, giver] <- give[j, giver] - result[["moved"]]
          }
        }
      }
      # rounding can leave values of about -1e-18
      transitions[j, , ] <- pmax(m, 0)
    }
    attr(transitions, "take") <- take
    attr(transitions, "give") <- give
    return(transitions)
  }

  # Moves target inflow from takers to givers (amounts per unit and type) in the grid cells where the
  # taker supplies the target and the giver has area left, in proportion to the area that can move.
  # Returns the inflow with the area that could not be moved per unit as attribute "remaining".
  .swapSources <- function(inflow, room, take, give, lrIndex) {
    remaining <- rep(0, nrow(take))
    for (j in which(rowSums(take) > 0)) {
      cells <- which(lrIndex == j)
      for (taker in which(take[j, ] > 0)) {
        for (giver in which(give[j, ] > 0)) {
          movable <- pmin(inflow[cells, taker], pmax(room[cells, giver], 0))
          amount <- min(take[j, taker], give[j, giver], sum(movable))
          if (amount > 0) {
            moved <- amount * movable / sum(movable)
            inflow[cells, taker] <- inflow[cells, taker] - moved
            inflow[cells, giver] <- inflow[cells, giver] + moved
            room[cells, giver] <- room[cells, giver] - moved
            take[j, taker] <- take[j, taker] - amount
            give[j, giver] <- give[j, giver] - amount
          }
        }
      }
      remaining[j] <- sum(take[j, ])
    }
    attr(inflow, "remaining") <- remaining
    return(inflow)
  }

  # One compensation step for a taker and a giver in one unit, moving at most amount; each step
  # is shared in proportion to the existing transitions:
  # 1. transitions from the taker to the giver are reduced;
  # 2. transitions from the taker to a third type are replaced by transitions from the giver;
  # 3. transitions from the giver to the taker are added, where allowed;
  # 4. transitions from the giver to a third type and from there to the taker are added, where
  #    allowed, in proportion to the area of the third types not yet used by outflows.
  .compensate <- function(step, m, taker, giver, amount, allowed, area) {
    if (step == 1) {
      moved <- min(amount, m[taker, giver])
      m[taker, giver] <- m[taker, giver] - moved
    } else if (step == 2) {
      third <- setdiff(which(allowed[giver, ] & m[taker, ] > 0), c(giver, taker))
      room <- sum(m[taker, third])
      moved <- min(amount, room)
      if (moved > 0) {
        share <- m[taker, third] / room
        m[taker, third] <- m[taker, third] - moved * share
        m[giver, third] <- m[giver, third] + moved * share
      }
    } else if (step == 3) {
      moved <- if (allowed[giver, taker]) amount else 0
      m[giver, taker] <- m[giver, taker] + moved
    } else {
      via <- setdiff(which(allowed[giver, ] & allowed[, taker]), c(giver, taker))
      unused <- pmax(area[via] - rowSums(m)[via], 0)
      moved <- min(amount, sum(unused))
      if (moved > 0) {
        share <- unused / sum(unused)
        m[giver, via] <- m[giver, via] + moved * share
        m[via, taker] <- m[via, taker] + moved * share
      }
    }
    return(list(m = m, moved = moved))
  }

  # New grid cell state: old state minus outflows plus inflows of the grid cell transitions ("from.to").
  .applyFlows <- function(state, flows) {
    from <- sub("[.].*$", "", colnames(flows))
    to <- sub("^.*[.]", "", colnames(flows))
    for (n in seq_len(ncol(flows))) {
      state[, from[n]] <- state[, from[n]] - flows[, n]
      state[, to[n]] <- state[, to[n]] + flows[, n]
    }

    # outflows can exceed the area by up to tol
    state[state < 0] <- 0
    return(state)
  }


  hrToLr <- as.character(mapping[["cluster"]][match(getCells(startHr), mapping[["cell"]])])
  .checkTransitionInput(startHr, hrToLr, statesLr, transitionsLr, constraints)
  types <- getNames(startHr)
  cells <- getCells(startHr)
  years <- getYears(transitionsLr)
  unitsLr <- getCells(statesLr)

  # unit states (units x types x years), unit transitions by year and stage (units x "from.to") and the
  # bounds of the constraints (grid cells x years)
  statesArray <- as.array(statesLr)[unitsLr, years, types, drop = FALSE]
  statesUnits <- aperm(statesArray, c(1, 3, 2))
  stageOf <- getItems(transitionsLr, dim = 3.1, full = TRUE)
  stagesLr <- list()
  for (year in years) {
    stagesLr[[year]] <- list()
    for (stage in unique(stageOf)) {
      items <- getNames(transitionsLr)[stageOf == stage]
      stagesLr[[year]][[stage]] <- matrix(transitionsLr[unitsLr, year, items], nrow = length(unitsLr),
                                          dimnames = list(unitsLr, sub("^[^.]*[.]", "", items)))
    }
  }
  bounds <- lapply(constraints, function(x) {
    bound <- if (is.null(x[["upper"]])) x[["lower"]] else x[["upper"]]
    return(matrix(bound[cells, , ], length(cells), dimnames = list(cells, getYears(bound))))
  })
  cellArea <- rowSums(.asMatrix(startHr))

  state <- .fitInitial(.asMatrix(startHr), hrToLr, .unitState(years[1], 0), tol)
  startOut <- .asMagpie(state, getYears(startHr))
  stateOut <- list()
  flowsOut <- list()
  diagnostics <- list()
  excess <- list()
  transitionChange <- list()
  previous <- NULL

  for (year in years) {
    # caps of the year from the constraints: the bound of the latest year at or before the year (the
    # first year before it, and the only one for bounds without year); a lower bound becomes an upper bound
    # on all other types
    caps <- list()
    for (name in names(constraints)) {
      x <- constraints[[name]]
      boundYears <- colnames(bounds[[name]])
      latest <- max(c(1, which(boundYears <= year)))
      bound <- bounds[[name]][, latest]
      if (is.null(x[["upper"]])) {
        caps[[name]] <- list(types = setdiff(types, x[["types"]]), cap = cellArea - bound)
      } else {
        caps[[name]] <- list(types = x[["types"]], cap = bound)
      }
      caps[[name]][["priority"]] <- if (is.null(x[["priority"]])) 1 else x[["priority"]]
    }

    flowsYear <- list()
    for (stage in names(stagesLr[[year]])) {
      input <- list()
      if (!is.null(stageRules)) {
        input <- as.list(stageRules(year, stage, .asMagpie(state, year), previous, layer))
      }
      for (weight in intersect(c("fromWeight", "toWeight"), names(input))) {
        stopifnot("weights need the grid cells of startHr in its order" = identical(getCells(input[[weight]]), cells))
        input[[weight]] <- .asMatrix(input[[weight]])
      }
      input[["caps"]] <- caps
      total <- stagesLr[[year]][[stage]]
      # a type with a target grid cell area gets its transitions from the target first; the other
      # transitions are balanced to the same unit states and allocated afterwards
      if (is.null(input[["targetArea"]])) {
        result <- .allocateStage(state, mapping, hrToLr, total, input, tol, maxIter)
      } else {
        if (ndata(input[["targetArea"]]) != 1) {
          stop("targetArea must name one type")
        }
        targetType <- getNames(input[["targetArea"]])
        targetHr <- as.vector(input[["targetArea"]][cells, , ])
        allowed <- NULL
        if (!is.null(input[["allowedTransitions"]])) {
          pairs <- input[["allowedTransitions"]]
          allowed <- matrix(FALSE, length(types), length(types), dimnames = list(types, types))
          allowed[cbind(sub("[.].*$", "", pairs), sub("^.*[.]", "", pairs))] <- TRUE
        }
        unitState <- .unitState(year, match(stage, names(stagesLr[[year]])))
        target <- .targetFlows(state, mapping, hrToLr, total, targetType, targetHr, unitState, allowed, tol)
        result <- .allocateStage(.applyFlows(state, target[["flows"]]), mapping, hrToLr, target[["rest"]],
                                 input, tol, maxIter)
        result[["flows"]] <- cbind(target[["flows"]], result[["flows"]])
        result[["transitionChange"]] <- target[["change"]]
      }
      state <- result[["state"]]
      flows <- result[["flows"]]
      colnames(flows) <- paste(rep(stage, ncol(flows)), colnames(flows), sep = ".")
      flowsYear[[stage]] <- flows
      diagnostics[[length(diagnostics) + 1]] <- data.frame(
        year = year, stage = stage, iterations = result[["iterations"]],
        maxResidual = max(c(0, result[["residual"]])), nExceeded = length(unique(result[["excess"]][["cluster"]]))
      )
      if (NROW(result[["excess"]]) > 0) {
        excess[[length(excess) + 1]] <- cbind(year = year, stage = stage, result[["excess"]])
      }
      if (NROW(result[["transitionChange"]]) > 0) {
        transitionChange[[length(transitionChange) + 1]] <- cbind(year = year, stage = stage,
                                                                  result[["transitionChange"]])
      }
    }
    stateOut[[year]] <- .asMagpie(state, year)
    # the grid cell transitions of the year, without those that are zero in all grid cells
    flowsYear <- do.call(cbind, unname(flowsYear))
    flowsOut[[year]] <- flowsYear[, colSums(flowsYear) > 0, drop = FALSE]
    previous <- NULL
    if (ncol(flowsOut[[year]]) > 0) {
      previous <- .asMagpie(flowsOut[[year]], year)
    }
    if (!is.null(updateLayer)) {
      layer <- updateLayer(year, .asMagpie(state, year), previous, layer)
    }
  }

  # grid cell states and transitions of all years, joined once at the end: each assignment into a magclass
  # object copies the whole object, so filling one object year by year takes several times its size in memory
  stateOut <- mbind(stateOut)

  # the transitions, the largest output, in the order of the stages, assembled in an array: assembled in magclass
  # (year by year or with mbind), they take five to seven times their size in memory
  items <- unique(unlist(lapply(flowsOut, colnames)))
  itemStage <- sub("[.].*$", "", items)
  items <- items[order(match(itemStage, unique(stageOf)))]
  values <- array(0, c(length(cells), length(years), length(items)), list(cells, years, items))
  for (year in years) {
    values[, year, colnames(flowsOut[[year]])] <- flowsOut[[year]]
    flowsOut[[year]] <- NULL
  }
  transitionsOut <- as.magpie(values, spatial = 1, temporal = 2)
  rm(values)

  # largest deviation of the unit sums of the grid cell states from the unit states
  deviation <- 0
  for (year in years) {
    sumLr <- rowsum(.asMatrix(stateOut[, year, ]), hrToLr)[unitsLr, , drop = FALSE]
    deviation <- max(deviation, abs(sumLr - statesUnits[, colnames(sumLr), year]))
  }

  return(list(startHr = startOut, stateHr = stateOut, transitionsHr = transitionsOut,
              diagnostics = do.call(rbind, diagnostics), excess = do.call(rbind, excess),
              transitionChange = do.call(rbind, transitionChange), layer = layer, deviationLr = deviation))
}
