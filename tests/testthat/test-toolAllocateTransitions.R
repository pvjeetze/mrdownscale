# helpers ------------------------------------------------------------------------------------

landTypes <- c("crop", "past", "forest", "urban", "other")

# synthetic land system: 3 clusters with 4 cells each; start = cluster sums of the initial state
makeLand <- function(seed = 1) {
  set.seed(seed)
  cluster <- rep(c("A", "B", "C"), each = 4)
  area <- runif(12, 0.5, 1)
  share <- matrix(runif(12 * 5, 0.1, 1), 12, 5)
  initial <- share / rowSums(share) * area
  dimnames(initial) <- list(paste0("c", 1:12), landTypes)
  list(initial = initial, cluster = cluster, start = rowsum(initial, cluster))
}

# cluster transitions: fraction of the source area in the start state, scaled by cluster
makeTotals <- function(start, fractions, scale = c(1, 0.8, 1.2)) {
  from <- sub("\\..*", "", names(fractions))
  out <- vapply(seq_along(fractions), function(n) start[, from[n]] * fractions[[n]] * scale,
                numeric(nrow(start)))
  dimnames(out) <- list(rownames(start), names(fractions))
  return(out)
}

# net change of the cluster states by transitions named "from.to"
netChange <- function(total, types = landTypes) {
  out <- matrix(0, nrow(total), length(types), dimnames = list(rownames(total), types))
  for (n in seq_len(ncol(total))) {
    ends <- strsplit(colnames(total)[n], ".", fixed = TRUE)[[1]]
    out[, ends[1]] <- out[, ends[1]] - total[, n]
    out[, ends[2]] <- out[, ends[2]] + total[, n]
  }
  return(out)
}

maxAbs <- function(a, b) max(abs(a - b))

# unit states at the end of each year (units x types x years) from the start and the transitions
unitStates <- function(start, transitions) {
  types <- colnames(start)
  out <- array(0, c(nrow(start), length(types), length(transitions)),
               list(rownames(start), types, names(transitions)))
  state <- start
  for (year in names(transitions)) {
    for (total in transitions[[year]]) state <- state + netChange(total, types)
    out[, , year] <- state
  }
  return(out)
}

# two years with two stages each
makeTransitions <- function(start) {
  year <- function(scale) {
    list(shift = makeTotals(start, c(forest.other = 0.1), scale),
         opt = makeTotals(start, c(crop.past = 0.1, past.forest = 0.08, crop.urban = 0.05,
                                   other.crop = 0.1, forest.crop = 0.05), scale))
  }
  list(y2000 = year(c(1, 0.8, 1.2)), y2010 = year(c(1.2, 1, 0.8)))
}

# stageRules that returns the same input in every stage
constantInput <- function(input) {
  function(year, stage, state, previous, layer) input
}

# stageRules with an urban target for stage "opt": start urban plus a gain that follows the past area
# and sums to the net urban change of the cluster
urbanTarget <- function(totals, cluster) {
  function(year, stage, state, previous, layer) {
    if (stage != "opt") {
      return(NULL)
    }
    weight <- state[, "past"]^2
    netUrban <- totals[, "crop.urban"] + totals[, "past.urban"]
    gain <- unname(netUrban[cluster]) * weight / rowsum(weight, cluster)[cluster, 1]
    list(targetArea = list(urban = state[, "urban"] + gain))
  }
}

# totals with the three transitions used for the target tests
makeTargetTotals <- function(start) {
  makeTotals(start, c(crop.urban = 0.05, past.urban = 0.03, crop.past = 0.1))
}

# tests --------------------------------------------------------------------------------------

test_that("cluster states and transitions are exact", {
  land <- makeLand()
  transitions <- makeTransitions(land[["start"]])
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                 unitStates(land[["start"]], transitions), transitions)

  expect_equal(dimnames(res[["stateHr"]])[[3]], c("start", "y2000", "y2010"))
  expect_true(all(res[["stateHr"]] >= 0))
  expect_lt(max(res[["diagnostics"]][["maxResidual"]]), 1e-8)

  expected <- land[["start"]]
  for (year in names(transitions)) {
    for (stage in names(transitions[[year]])) {
      total <- transitions[[year]][[stage]]
      flows <- res[["transitionsHr"]][[year]][[stage]]
      expected <- expected + netChange(total)
      expect_setequal(colnames(flows), colnames(total))
      expect_lt(maxAbs(rowsum(flows, land[["cluster"]])[, colnames(flows)], total[, colnames(flows)]), 1e-8)
    }
    expect_lt(maxAbs(rowsum(res[["stateHr"]][, , year], land[["cluster"]]), expected), 1e-8)
    expect_lt(maxAbs(rowSums(res[["stateHr"]][, , year]), rowSums(land[["initial"]])), 1e-8)
  }
})

test_that("no cell gains more than the inflow of its cluster", {
  land <- makeLand()
  transitions <- makeTransitions(land[["start"]])
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                 unitStates(land[["start"]], transitions), transitions)

  previous <- "start"
  for (year in names(transitions)) {
    gain <- rowsum(pmax(res[["stateHr"]][, , year] - res[["stateHr"]][, , previous], 0), land[["cluster"]])
    inflow <- matrix(0, 3, length(landTypes), dimnames = dimnames(gain))
    for (total in transitions[[year]]) {
      for (n in seq_len(ncol(total))) {
        to <- strsplit(colnames(total)[n], ".", fixed = TRUE)[[1]][2]
        inflow[, to] <- inflow[, to] + total[, n]
      }
    }
    expect_lte(max(gain - inflow), 1e-8)
    # the bound is not trivial: the cells gain a positive area
    expect_gt(min(gain[, "past"]), 0)
    previous <- year
  }
})

test_that("transitions out of a type are proportional to its area by default", {
  land <- makeLand()
  total <- makeTotals(land[["start"]], c(crop.other = 0.2))
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                 unitStates(land[["start"]], list(y2000 = list(opt = total))),
                                 list(y2000 = list(opt = total)))

  flow <- res[["transitionsHr"]][["y2000"]][["opt"]][, "crop.other"]
  cropSum <- rowsum(land[["initial"]][, "crop"], land[["cluster"]])[land[["cluster"]], 1]
  expected <- total[land[["cluster"]], "crop.other"] * land[["initial"]][, "crop"] / cropSum
  expect_equal(unname(flow), unname(expected), tolerance = 1e-10)
})

test_that("initial state is fitted to start", {
  land <- makeLand()
  # perturb the types within each cluster so that the cluster areas stay the same
  shift <- c(1, -1, 1, -1, 0) * 1e-3
  start <- land[["start"]] + rbind(shift, -shift, shift)
  total <- makeTotals(start, c(crop.other = 0.2))
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                 unitStates(start, list(y2000 = list(opt = total))),
                                 list(y2000 = list(opt = total)))

  expect_lt(maxAbs(rowsum(res[["stateHr"]][, , "start"], land[["cluster"]]), start), 1e-10)
  expect_lt(maxAbs(rowSums(res[["stateHr"]][, , "start"]), rowSums(land[["initial"]])), 1e-8)
  expect_true(all(res[["stateHr"]] >= 0))
  # the fit changes the initial state
  expect_gt(maxAbs(res[["stateHr"]][, , "start"], land[["initial"]]), 1e-5)
})

test_that("weights from stageRules change the placement", {
  land <- makeLand()
  transitions <- list(y2000 = list(opt = makeTotals(land[["start"]], c(forest.other = 0.1))))
  plain <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                   unitStates(land[["start"]], transitions), transitions)

  weight <- matrix(0.01, 12, 1, dimnames = list(NULL, "other"))
  first <- c(1, 5, 9)
  weight[first, ] <- 1
  weighted <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                      unitStates(land[["start"]], transitions), transitions,
                                      stageRules = constantInput(list(toWeight = weight)))

  flowPlain <- plain[["transitionsHr"]][["y2000"]][["opt"]][, "forest.other"]
  flowWeighted <- weighted[["transitionsHr"]][["y2000"]][["opt"]][, "forest.other"]
  expect_gt(sum(flowWeighted[first]), sum(flowPlain[first]) + 0.1)
  total <- transitions[["y2000"]][["opt"]][, "forest.other"]
  expect_lt(maxAbs(rowsum(flowWeighted, land[["cluster"]])[, 1], total), 1e-8)
})

test_that("soft caps hold where feasible and are exceeded by the shortfall where not", {
  land <- makeLand()
  total <- makeTotals(land[["start"]], c(forest.other = 0.3))
  demand <- total[, "forest.other"]
  # A: feasible and binding in some cells, B: caps sum to half the demand, C: no cap
  cap <- c(demand[["A"]] * c(0.1, 0.2, 0.3, 0.5), rep(demand[["B"]] / 8, 4), rep(Inf, 4))
  constraints <- list(intoOther = list(types = "other", upper = land[["initial"]][, "other"] + cap))
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                 unitStates(land[["start"]], list(y2000 = list(opt = total))),
                                 list(y2000 = list(opt = total)), constraints = constraints)

  flow <- res[["transitionsHr"]][["y2000"]][["opt"]][, "forest.other"]
  cellsA <- land[["cluster"]] == "A"
  expect_lte(max(flow[cellsA] - cap[cellsA]), 1e-6)
  expect_lt(min(abs(flow[cellsA] - cap[cellsA])), 1e-6)
  # in B the caps are exceeded by the shortfall only, in A by less than 1e-5 of the flow
  excess <- setNames(res[["excess"]][["area"]], res[["excess"]][["cluster"]])
  expect_equal(excess[["B"]], demand[["B"]] / 2, tolerance = 1e-4)
  expect_lt(sum(excess[names(excess) == "A"]), 1e-5 * demand[["A"]])
  expect_false("C" %in% names(excess))
  expect_lt(maxAbs(rowsum(flow, land[["cluster"]])[, 1], demand), 1e-8)
  expect_lt(maxAbs(rowsum(res[["stateHr"]][, , "y2000"], land[["cluster"]]),
                   land[["start"]] + netChange(total)), 1e-8)
})

test_that("cells above a cap give area first", {
  land <- makeLand()
  total <- makeTotals(land[["start"]], c(crop.other = 0.3))
  demand <- total[land[["cluster"]], "crop.other"]
  # the first cell of each cluster is above its crop cap by 60 % of the cluster demand
  first <- c(1, 5, 9)
  above <- rep(0, 12)
  above[first] <- pmin(0.6 * demand[first], land[["initial"]][first, "crop"])
  constraints <- list(crop = list(types = "crop", upper = land[["initial"]][, "crop"] - above))
  plain <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                   unitStates(land[["start"]], list(y2000 = list(opt = total))),
                                   list(y2000 = list(opt = total)))
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                 unitStates(land[["start"]], list(y2000 = list(opt = total))),
                                 list(y2000 = list(opt = total)),
                                 constraints = constraints)

  flow <- res[["transitionsHr"]][["y2000"]][["opt"]][, "crop.other"]
  flowPlain <- plain[["transitionsHr"]][["y2000"]][["opt"]][, "crop.other"]
  expect_true(all(flowPlain[first] < above[first]))
  expect_true(all(flow[first] >= above[first] - 1e-6))
  expect_lt(maxAbs(rowsum(flow, land[["cluster"]])[, 1], total[, "crop.other"]), 1e-8)
  # an outflow does not add to the area above the cap
  expect_null(res[["excess"]])
})

test_that("caps of lower priority are exceeded first", {
  land <- makeLand()
  total <- makeTotals(land[["start"]], c(forest.other = 0.3))
  # in each cluster, the first two cells are closed by a cap of low priority and the last two by a
  # cap of high priority, so one of the caps has to be exceeded
  # (forest kept at its area by a lower bound)
  closed <- rep(c(TRUE, TRUE, FALSE, FALSE), 3)
  forestArea <- land[["initial"]][, "forest"]
  constraints <- list(
    low = list(types = "forest", lower = ifelse(closed, forestArea, 0), priority = 1),
    high = list(types = "forest", lower = ifelse(closed, 0, forestArea), priority = 2)
  )
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                 unitStates(land[["start"]], list(y2000 = list(opt = total))),
                                 list(y2000 = list(opt = total)), constraints = constraints)

  flow <- res[["transitionsHr"]][["y2000"]][["opt"]][, "forest.other"]
  expect_lt(sum(flow[!closed]), 1e-4 * sum(flow))
  excess <- tapply(res[["excess"]][["area"]], res[["excess"]][["constraint"]], sum)
  expect_equal(excess[["low"]], sum(flow), tolerance = 1e-4)
  expect_lt(sum(excess[names(excess) == "high"]), 1e-4 * sum(flow))
  expect_lt(maxAbs(rowsum(flow, land[["cluster"]])[, 1], total[, "forest.other"]), 1e-8)
})

test_that("a layer is updated at the end of each year and passed to the next year", {
  land <- makeLand()
  transitions <- list(y2000 = list(opt = matrix(0.1, 3, 1, dimnames = list(c("A", "B", "C"), "crop.forest"))),
                      y2010 = list(opt = matrix(0.2, 3, 1, dimnames = list(c("A", "B", "C"), "crop.forest"))))
  # the layer adds up the cell flows of each year; stageRules records the layer it receives
  seen <- new.env()
  input <- function(year, stage, state, previous, layer) {
    assign(year, layer, envir = seen)
    return(NULL)
  }
  updateLayer <- function(year, state, flows, layer) layer + rowSums(flows[["opt"]])
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                 unitStates(land[["start"]], transitions), transitions,
                                 stageRules = input, layer = rep(0, nrow(land[["initial"]])),
                                 updateLayer = updateLayer)
  expect_equal(seen[["y2000"]], rep(0, nrow(land[["initial"]])))
  expect_equal(unname(seen[["y2010"]]), unname(rowSums(res[["transitionsHr"]][["y2000"]][["opt"]])))
  expect_equal(unname(res[["layer"]]), unname(rowSums(res[["transitionsHr"]][["y2000"]][["opt"]]) +
                                                rowSums(res[["transitionsHr"]][["y2010"]][["opt"]])))
})

test_that("a target type follows its cell target and keeps the cluster states", {
  land <- makeLand()
  total <- makeTargetTotals(land[["start"]])
  input <- urbanTarget(total, land[["cluster"]])
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                 unitStates(land[["start"]], list(y2000 = list(opt = total))),
                                 list(y2000 = list(opt = total)), stageRules = input)

  target <- input("y2000", "opt", land[["initial"]], NULL)[["targetArea"]][["urban"]]
  expect_lt(maxAbs(res[["stateHr"]][, "urban", "y2000"], target), 1e-8)
  expect_lt(maxAbs(rowsum(res[["stateHr"]][, , "y2000"], land[["cluster"]]),
                   land[["start"]] + netChange(total)), 1e-8)
  expect_true(all(res[["stateHr"]] >= 0))
  expect_gt(nrow(res[["transitionChange"]]), 0)
  expect_lt(maxAbs(rowSums(res[["stateHr"]][, , "y2000"]), rowSums(land[["initial"]])), 1e-8)
})

# state at the end of one stage "opt" with a target area for one type (one unit A)
targetStage <- function(state, total, targetArea, allowed) {
  state <- withCells(state)
  startLr <- matrix(colSums(state), 1, dimnames = list("A", colnames(state)))
  rules <- function(year, stage, stateHr, previous, layer) {
    return(list(targetArea = targetArea, allowedTransitions = allowed))
  }
  res <- toolAllocateTransitions(state, cellMapping(rep("A", nrow(state))),
                                 unitStates(startLr, list(y2000 = list(opt = total))),
                                 list(y2000 = list(opt = total)), rules)
  return(res[["stateHr"]][, , "y2000"])
}

test_that("target inflow moves back to the giving type where compensation is not possible", {
  # cell 1 has forest and other land, cell 2 only other land; the model takes urban land from
  # forest, the target places it in cell 1, from forest and other land in proportion to their area;
  # forest cannot go to other land, so the target takes forest instead of other land in cell 1
  types <- c("forest", "other", "urban")
  state <- matrix(c(0.5, 0.5, 0, 0, 1, 0), 2, byrow = TRUE, dimnames = list(NULL, types))
  total <- matrix(0.1, 1, 1, dimnames = list("A", "forest.urban"))
  allowed <- matrix(FALSE, 3, 3, dimnames = list(types, types))
  allowed[c("forest", "other"), "urban"] <- TRUE
  allowed["other", "forest"] <- TRUE
  res <- targetStage(state, total, list(urban = c(0.1, 0)), allowed)
  expect_equal(unname(res[, "urban"]), c(0.1, 0), tolerance = 1e-10)
  expect_equal(colSums(res), c(forest = 0.4, other = 1.5, urban = 0.1), tolerance = 1e-10)
})

test_that("a cell keeps the target type where no type in the cell can take it", {
  # urban shrinks in cell 1 by the target, but cell 1 holds only urban land and forest, which cannot
  # gain; the model has no outflow of urban land
  types <- c("forest", "other", "urban")
  state <- matrix(c(0.5, 0, 0.5, 0, 1, 0), 2, byrow = TRUE, dimnames = list(NULL, types))
  total <- matrix(0.1, 1, 1, dimnames = list("A", "other.urban"))
  allowed <- matrix(FALSE, 3, 3, dimnames = list(types, types))
  allowed["other", "urban"] <- TRUE
  allowed["forest", "other"] <- TRUE
  res <- targetStage(state, total, list(urban = c(0.3, 0.1)), allowed)
  expect_equal(colSums(res), c(forest = 0.5, other = 0.9, urban = 0.6), tolerance = 1e-10)
  expect_equal(unname(res[, "urban"]), c(0.5, 0.1), tolerance = 1e-10)

  # with a tiny area of a type that can gain, the cell loses the target type to it, and the cluster
  # states still hold
  state[1, "other"] <- 1e-33
  res <- expect_silent(targetStage(state, total, list(urban = c(0.3, 0.1)), allowed))
  expect_equal(colSums(res), c(forest = 0.5, other = 0.9, urban = 0.6), tolerance = 1e-10)
})

test_that("compensation respects allowed", {
  land <- makeLand()
  total <- makeTargetTotals(land[["start"]])
  target <- urbanTarget(total, land[["cluster"]])
  allowedBut <- function(forbidden) {
    allowed <- matrix(TRUE, 5, 5, dimnames = list(landTypes, landTypes))
    allowed[forbidden] <- FALSE
    return(allowed)
  }
  run <- function(allowed) {
    input <- function(year, stage, state, previous, layer) {
      out <- target(year, stage, state, previous, layer)
      if (!is.null(out)) {
        out[["allowedTransitions"]] <- allowed
      }
      return(out)
    }
    toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                            unitStates(land[["start"]], list(y2000 = list(opt = total))),
                            list(y2000 = list(opt = total)),
                            stageRules = input)
  }
  area <- function(res, transition) {
    flows <- res[["transitionsHr"]][["y2000"]][["opt"]]
    if (transition %in% colnames(flows)) sum(flows[, transition]) else 0
  }
  expected <- land[["start"]] + netChange(total)

  # crop.forest is added by the balancing if it is allowed
  expect_gt(area(run(allowedBut(NULL)), "crop.forest"), 1e-3)
  res <- run(allowedBut(cbind("crop", "forest")))
  expect_equal(area(res, "crop.forest"), 0)
  expect_lt(maxAbs(rowsum(res[["stateHr"]][, , "y2000"], land[["cluster"]]), expected), 1e-8)
  expect_gt(nrow(res[["transitionChange"]]), 0)

  # without crop.forest, crop.other and past.other the other transitions cannot compensate; the
  # target then takes the giving types instead, and the cluster states hold
  forbidden <- rbind(c("crop", "forest"), c("crop", "other"), c("past", "other"))
  expect_silent(res <- run(allowedBut(forbidden)))
  for (transition in c("crop.forest", "crop.other", "past.other")) {
    expect_equal(area(res, transition), 0)
  }
  expect_lt(maxAbs(rowsum(res[["stateHr"]][, , "y2000"], land[["cluster"]]), expected), 1e-8)
})

test_that("input errors", {
  land <- makeLand()
  total <- makeTotals(land[["start"]], c(crop.other = 0.2))
  good <- list(y2000 = list(opt = total))
  states <- unitStates(land[["start"]], good)
  run <- function(initial = land[["initial"]], statesLr = states, stages = good) {
    toolAllocateTransitions(initial, cellMapping(land[["cluster"]]), statesLr, stages)
  }
  expect_silent(run())

  bad <- land[["initial"]]
  colnames(bad)[1] <- "crop.land"
  expect_error(run(initial = bad), "must not contain '[.]'")

  bad <- land[["initial"]]
  bad[1, 1] <- -0.1
  expect_error(run(initial = bad), "startHr must not contain NA or negative values")

  expect_error(run(statesLr = states[, -1, , drop = FALSE]), "statesLr must be an array")
  expect_error(run(statesLr = states[, , 1, drop = FALSE], stages = list(y2010 = good[[1]])),
               "statesLr must be an array")

  bad <- total
  colnames(bad) <- "crop.water"
  expect_error(run(stages = list(y2000 = list(opt = bad))), "with types of startHr")

  bad <- total
  rownames(bad) <- c("A", "B", "X")
  expect_error(run(stages = list(y2000 = list(opt = bad))), "must have the units of statesLr as row names")

  runWith <- function(constraint) {
    toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                            unitStates(land[["start"]], good), good,
                            constraints = list(a = constraint))
  }
  expect_error(runWith(list(types = "water", upper = rep(1, 12))), "contains types not in startHr")
  expect_error(runWith(list(types = "crop", upper = rep(1, 12), lower = rep(0, 12))), "exactly one of upper and lower")
  expect_error(runWith(list(types = "crop", upper = rep(1, 5))), "one value per grid cell")
  expect_error(runWith(list(types = "crop", upper = matrix(1, 12, 2))), "one value per grid cell")
})

test_that("constraints use the bound of the latest year, and a lower bound limits the other types", {
  startHr <- matrix(c(0, 10, 0, 10), 2, byrow = TRUE, dimnames = list(c("c1", "c2"), c("crop", "other")))
  startLr <- matrix(c(0, 20), 1, dimnames = list("A", c("crop", "other")))
  expansion <- list(opt = matrix(2, 1, 1, dimnames = list("A", "other.crop")))
  # cropland bound by year: y2015 uses y2010, years before y2000 use y2000
  bound <- matrix(c(10, 10, 0.5, 10, 10, 0.5), 2, dimnames = list(NULL, c("y2000", "y2010", "y2020")))
  upper <- list(crop = list(types = "crop", upper = bound))
  res <- toolAllocateTransitions(startHr, cellMapping(c("A", "A")),
                                 unitStates(startLr, list(y2015 = expansion)), list(y2015 = expansion),
                                 constraints = upper)
  expect_equal(unname(res[["stateHr"]][, "crop", "y2015"]), c(0.5, 1.5), tolerance = 1e-6)
  res <- toolAllocateTransitions(startHr, cellMapping(c("A", "A")),
                                 unitStates(startLr, list(y1995 = expansion)), list(y1995 = expansion),
                                 constraints = upper)
  expect_equal(unname(res[["stateHr"]][, "crop", "y1995"]), c(1, 1), tolerance = 1e-6)
  # a lower bound on other land keeps the cropland of cell 1 within the rest of the cell
  lower <- list(other = list(types = "other", lower = c(9.5, 0)))
  res <- toolAllocateTransitions(startHr, cellMapping(c("A", "A")),
                                 unitStates(startLr, list(y2015 = expansion)), list(y2015 = expansion),
                                 constraints = lower)
  expect_equal(unname(res[["stateHr"]][, "crop", "y2015"]), c(0.5, 1.5), tolerance = 1e-6)
})

test_that("the deviation from the unit states is reported, and a target is scaled to the unit totals", {
  land <- makeLand()
  transitions <- makeTransitions(land[["start"]])
  states <- unitStates(land[["start"]], transitions)
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]), states, transitions)
  expect_lt(res[["deviationLr"]], 1e-8)
  # unit states that differ from the transitions by 0.01 show as deviation
  shifted <- states
  shifted[, "crop", "y2010"] <- shifted[, "crop", "y2010"] + 0.01
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]), shifted, transitions)
  expect_equal(res[["deviationLr"]], 0.01, tolerance = 1e-6)

  # a target given as distribution (three times the area) gives the same result as the target in Mha
  total <- makeTargetTotals(land[["start"]])
  stages <- list(y2000 = list(opt = total))
  rules <- urbanTarget(total, land[["cluster"]])
  scaled <- function(year, stage, state, previous, layer) {
    input <- rules(year, stage, state, previous, layer)
    if (!is.null(input)) {
      input[["targetArea"]][["urban"]] <- 3 * input[["targetArea"]][["urban"]]
    }
    return(input)
  }
  inMha <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                   unitStates(land[["start"]], stages), stages, rules)
  asShare <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                     unitStates(land[["start"]], stages), stages, scaled)
  expect_lt(maxAbs(inMha[["stateHr"]], asShare[["stateHr"]]), 1e-8)
})

test_that("magclass bounds are matched by grid cell", {
  startHr <- matrix(c(0, 10, 0, 10), 2, byrow = TRUE, dimnames = list(c("c1", "c2"), c("crop", "other")))
  startLr <- matrix(c(0, 20), 1, dimnames = list("A", c("crop", "other")))
  stages <- list(y2015 = list(opt = matrix(2, 1, 1, dimnames = list("A", "other.crop"))))
  bound <- matrix(c(0.5, 10), 2, dimnames = list(NULL, "y2010"))
  asMatrix <- toolAllocateTransitions(startHr, cellMapping(c("A", "A")), unitStates(startLr, stages), stages,
                                      constraints = list(crop = list(types = "crop", upper = bound)))
  # the same bound as magclass object with the grid cells in reverse order
  boundMagpie <- magclass::new.magpie(c("c2", "c1"), "y2010", "crop", fill = c(10, 0.5))
  asMagpie <- toolAllocateTransitions(startHr, cellMapping(c("A", "A")), unitStates(startLr, stages), stages,
                                      constraints = list(crop = list(types = "crop", upper = boundMagpie)))
  expect_equal(asMagpie[["stateHr"]], asMatrix[["stateHr"]])
  expect_equal(unname(asMagpie[["stateHr"]][, "crop", "y2015"]), c(0.5, 1.5), tolerance = 1e-6)
})
