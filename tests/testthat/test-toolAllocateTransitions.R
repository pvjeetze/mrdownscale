# helpers ------------------------------------------------------------------------------------

landTypes <- c("crop", "past", "forest", "urban", "other")

# sums over the grid cells of each unit
unitSums <- function(x, cluster) {
  return(toolAggregate(x, cellMapping(cluster), from = "cell", to = "cluster"))
}

maxAbs <- function(a, b) {
  return(max(abs(as.vector(a) - as.vector(b))))
}

# synthetic land system: 3 clusters with 4 grid cells each; start = cluster sums of the initial state
makeLand <- function(seed = 1) {
  set.seed(seed)
  cluster <- rep(c("A", "B", "C"), each = 4)
  area <- runif(12, 0.5, 1)
  share <- matrix(runif(12 * 5, 0.1, 1), 12, 5)
  initial <- new.magpie(paste0("c", 1:12), "y1995", landTypes, fill = share / rowSums(share) * area)
  list(initial = initial, cluster = cluster, start = unitSums(initial, cluster))
}

# unit totals of transitions "from.to": fraction of the source area in the start state, scaled by unit
makeTotals <- function(start, fractions, scale = c(1, 0.8, 1.2)) {
  from <- sub("[.].*$", "", names(fractions))
  out <- new.magpie(getCells(start), NULL, names(fractions), fill = 0)
  for (n in seq_along(fractions)) {
    out[, , n] <- setYears(start[, , from[n]], NULL) * fractions[[n]] * scale
  }
  return(out)
}

# unit transitions of one stage "opt" in one year
oneStage <- function(total, year = "y2000") {
  return(add_dimension(setYears(total, year), dim = 3.1, add = "stage", nm = "opt"))
}

# two years with two stages each
makeTransitions <- function(start) {
  .year <- function(year, scale) {
    shift <- add_dimension(makeTotals(start, c(forest.other = 0.1), scale), dim = 3.1, add = "stage", nm = "shift")
    opt <- add_dimension(makeTotals(start, c(crop.past = 0.1, past.forest = 0.08, crop.urban = 0.05,
                                             other.crop = 0.1, forest.crop = 0.05), scale),
                         dim = 3.1, add = "stage", nm = "opt")
    return(setYears(mbind(shift, opt), year))
  }
  return(mbind(.year("y2000", c(1, 0.8, 1.2)), .year("y2010", c(1.2, 1, 0.8))))
}

# net change of the unit states by the transitions of one year ("stage.from.to" or "from.to")
netChange <- function(transitions, types = landTypes) {
  out <- new.magpie(getCells(transitions), NULL, types, fill = 0)
  pairs <- sub("^[^.]+[.]([^.]+[.][^.]+)$", "\\1", getNames(transitions))
  for (n in seq_along(pairs)) {
    flow <- collapseDim(setYears(transitions[, , n], NULL), dim = 3)
    out[, , sub("[.].*$", "", pairs[n])] <- out[, , sub("[.].*$", "", pairs[n])] - flow
    out[, , sub("^.*[.]", "", pairs[n])] <- out[, , sub("^.*[.]", "", pairs[n])] + flow
  }
  return(out)
}

# unit states at the end of each year from the start and the transitions
unitStates <- function(start, transitions) {
  out <- new.magpie(getCells(start), getYears(transitions), getNames(start), fill = 0)
  state <- setYears(start, NULL)
  for (year in getYears(transitions)) {
    state <- state + netChange(transitions[, year, ], getNames(start))
    out[, year, ] <- state
  }
  return(out)
}

# stageRules that returns the same input in every stage
constantInput <- function(input) {
  function(year, stage, state, previous, layer) input
}

# stageRules with an urban target for stage "opt": start urban plus a gain that follows the past area
# and sums to the net urban change of the unit
urbanTarget <- function(totals, cluster) {
  function(year, stage, state, previous, layer) {
    if (stage != "opt") {
      return(NULL)
    }
    weight <- as.vector(state[, , "past"])^2
    netUrban <- as.vector(dimSums(totals[, , c("crop.urban", "past.urban")], dim = 3))
    gain <- netUrban[match(cluster, getCells(totals))] * weight / as.vector(tapply(weight, cluster, sum)[cluster])
    list(targetArea = state[, , "urban"] + gain)
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
  states <- unitStates(land[["start"]], transitions)
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]), states, transitions)

  expect_equal(getYears(res[["stateHr"]]), c("y2000", "y2010"))
  expect_true(all(res[["stateHr"]] >= 0))
  expect_lt(max(res[["diagnostics"]][["maxResidual"]]), 1e-8)
  expect_setequal(getNames(res[["transitionsHr"]]), getNames(transitions))
  flowsLr <- unitSums(res[["transitionsHr"]], land[["cluster"]])
  expect_lt(maxAbs(flowsLr[getCells(transitions), , getNames(transitions)], transitions), 1e-8)
  expect_lt(maxAbs(unitSums(res[["stateHr"]], land[["cluster"]]), states), 1e-8)
  # the area of each grid cell stays the same
  area <- as.vector(dimSums(land[["initial"]], dim = 3))
  expect_lt(maxAbs(dimSums(res[["stateHr"]], dim = 3), rep(area, 2)), 1e-8)
})

test_that("no cell gains more than the inflow of its cluster", {
  land <- makeLand()
  transitions <- makeTransitions(land[["start"]])
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                 unitStates(land[["start"]], transitions), transitions)

  to <- sub("^.*[.]", "", getNames(transitions))
  previous <- setYears(res[["startHr"]], NULL)
  for (year in getYears(transitions)) {
    current <- setYears(res[["stateHr"]][, year, ], NULL)
    gain <- unitSums(pmax(current - previous, 0), land[["cluster"]])
    inflow <- gain * 0
    for (k in unique(to)) {
      inflow[, , k] <- dimSums(transitions[, year, to == k], dim = 3)
    }
    expect_lte(max(gain - inflow), 1e-8)
    # the bound is not trivial: the cells gain a positive area
    expect_gt(min(gain[, , "past"]), 0)
    previous <- current
  }
})

test_that("transitions out of a type are proportional to its area by default", {
  land <- makeLand()
  total <- makeTotals(land[["start"]], c(crop.other = 0.2))
  transitions <- oneStage(total)
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                 unitStates(land[["start"]], transitions), transitions)

  crop <- as.vector(land[["initial"]][, , "crop"])
  expected <- as.vector(total)[match(land[["cluster"]], getCells(total))] * crop /
    as.vector(tapply(crop, land[["cluster"]], sum)[land[["cluster"]]])
  expect_equal(as.vector(res[["transitionsHr"]][, , "opt.crop.other"]), expected, tolerance = 1e-10)
})

test_that("initial state is fitted to start", {
  land <- makeLand()
  # perturb the types within each cluster so that the cluster areas stay the same
  shift <- c(1, -1, 1, -1, 0) * 1e-3
  start <- land[["start"]]
  start["A", , ] <- start["A", , ] + shift
  start["B", , ] <- start["B", , ] - shift
  start["C", , ] <- start["C", , ] + shift
  transitions <- oneStage(makeTotals(start, c(crop.other = 0.2)))
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                 unitStates(start, transitions), transitions)

  expect_lt(maxAbs(unitSums(res[["startHr"]], land[["cluster"]]), start), 1e-10)
  expect_lt(maxAbs(dimSums(res[["startHr"]], dim = 3), dimSums(land[["initial"]], dim = 3)), 1e-8)
  expect_true(all(res[["stateHr"]] >= 0))
  # the fit changes the initial state
  expect_gt(maxAbs(res[["startHr"]], land[["initial"]]), 1e-5)
})

test_that("weights from stageRules change the placement", {
  land <- makeLand()
  transitions <- oneStage(makeTotals(land[["start"]], c(forest.other = 0.1)))
  states <- unitStates(land[["start"]], transitions)
  plain <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]), states, transitions)

  first <- c("c1", "c5", "c9")
  weight <- new.magpie(getCells(land[["initial"]]), NULL, "other", fill = 0.01)
  weight[first, , ] <- 1
  weighted <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]), states, transitions,
                                      stageRules = constantInput(list(toWeight = weight)))

  flowPlain <- plain[["transitionsHr"]][first, , "opt.forest.other"]
  flowWeighted <- weighted[["transitionsHr"]][first, , "opt.forest.other"]
  expect_gt(sum(flowWeighted), sum(flowPlain) + 0.1)
  expect_lt(maxAbs(unitSums(weighted[["transitionsHr"]], land[["cluster"]]), transitions), 1e-8)
})

test_that("soft caps hold where feasible and are exceeded by the shortfall where not", {
  land <- makeLand()
  total <- makeTotals(land[["start"]], c(forest.other = 0.3))
  transitions <- oneStage(total)
  demand <- stats::setNames(as.vector(total), getCells(total))
  # A: feasible and binding in some cells, B: caps sum to half the demand, C: no cap
  cap <- c(demand[["A"]] * c(0.1, 0.2, 0.3, 0.5), rep(demand[["B"]] / 8, 4), rep(Inf, 4))
  constraints <- list(intoOther = list(types = "other", upper = land[["initial"]][, , "other"] + cap))
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                 unitStates(land[["start"]], transitions), transitions, constraints = constraints)

  flow <- as.vector(res[["transitionsHr"]][, , "opt.forest.other"])
  cellsA <- land[["cluster"]] == "A"
  expect_lte(max(flow[cellsA] - cap[cellsA]), 1e-6)
  expect_lt(min(abs(flow[cellsA] - cap[cellsA])), 1e-6)
  # in B the caps are exceeded by the shortfall only, in A by less than 1e-5 of the flow
  excess <- stats::setNames(res[["excess"]][["area"]], res[["excess"]][["cluster"]])
  expect_equal(excess[["B"]], demand[["B"]] / 2, tolerance = 1e-4)
  expect_lt(sum(excess[names(excess) == "A"]), 1e-5 * demand[["A"]])
  expect_false("C" %in% names(excess))
  expect_lt(maxAbs(tapply(flow, land[["cluster"]], sum), demand), 1e-8)
  expect_lt(maxAbs(unitSums(res[["stateHr"]], land[["cluster"]]), unitStates(land[["start"]], transitions)), 1e-8)
})

test_that("cells above a cap give area first", {
  land <- makeLand()
  total <- makeTotals(land[["start"]], c(crop.other = 0.3))
  transitions <- oneStage(total)
  states <- unitStates(land[["start"]], transitions)
  demand <- as.vector(total)[match(land[["cluster"]], getCells(total))]
  crop <- as.vector(land[["initial"]][, , "crop"])
  # the first cell of each cluster is above its crop cap by 60 % of the cluster demand
  first <- c(1, 5, 9)
  above <- rep(0, 12)
  above[first] <- pmin(0.6 * demand[first], crop[first])
  constraints <- list(crop = list(types = "crop", upper = land[["initial"]][, , "crop"] - above))
  plain <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]), states, transitions)
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]), states, transitions,
                                 constraints = constraints)

  flow <- as.vector(res[["transitionsHr"]][, , "opt.crop.other"])
  flowPlain <- as.vector(plain[["transitionsHr"]][, , "opt.crop.other"])
  expect_true(all(flowPlain[first] < above[first]))
  expect_true(all(flow[first] >= above[first] - 1e-6))
  expect_lt(maxAbs(tapply(flow, land[["cluster"]], sum), total), 1e-8)
  # an outflow does not add to the area above the cap
  expect_null(res[["excess"]])
})

test_that("caps of lower priority are exceeded first", {
  land <- makeLand()
  transitions <- oneStage(makeTotals(land[["start"]], c(forest.other = 0.3)))
  # in each cluster, the first two cells are closed by a cap of low priority and the last two by a
  # cap of high priority, so one of the caps has to be exceeded
  # (forest kept at its area by a lower bound)
  closed <- rep(c(TRUE, TRUE, FALSE, FALSE), 3)
  forestArea <- land[["initial"]][, , "forest"]
  constraints <- list(low = list(types = "forest", lower = forestArea * closed, priority = 1),
                      high = list(types = "forest", lower = forestArea * !closed, priority = 2))
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                 unitStates(land[["start"]], transitions), transitions, constraints = constraints)

  flow <- as.vector(res[["transitionsHr"]][, , "opt.forest.other"])
  expect_lt(sum(flow[!closed]), 1e-4 * sum(flow))
  excess <- tapply(res[["excess"]][["area"]], res[["excess"]][["constraint"]], sum)
  expect_equal(excess[["low"]], sum(flow), tolerance = 1e-4)
  expect_lt(sum(excess[names(excess) == "high"]), 1e-4 * sum(flow))
  expect_lt(maxAbs(tapply(flow, land[["cluster"]], sum), transitions), 1e-8)
})

test_that("a layer is updated at the end of each year and passed to the next year", {
  land <- makeLand()
  transitions <- new.magpie(c("A", "B", "C"), c("y2000", "y2010"), "opt.crop.forest",
                            fill = c(0.1, 0.1, 0.1, 0.2, 0.2, 0.2))
  # the layer adds up the cell flows of each year; stageRules records the layer it receives
  seen <- new.env()
  input <- function(year, stage, state, previous, layer) {
    assign(year, layer, envir = seen)
    return(NULL)
  }
  updateLayer <- function(year, state, flows, layer) {
    return(layer + setYears(dimSums(flows, dim = 3), NULL))
  }
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                 unitStates(land[["start"]], transitions), transitions,
                                 stageRules = input, layer = 0, updateLayer = updateLayer)
  expect_equal(seen[["y2000"]], 0)
  flows <- res[["transitionsHr"]]
  expect_equal(as.vector(seen[["y2010"]]), as.vector(dimSums(flows[, "y2000", ], dim = 3)))
  expect_equal(as.vector(res[["layer"]]), as.vector(dimSums(flows, dim = c(2, 3))))
})

test_that("a target type follows its cell target and keeps the cluster states", {
  land <- makeLand()
  total <- makeTargetTotals(land[["start"]])
  transitions <- oneStage(total)
  input <- urbanTarget(total, land[["cluster"]])
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                 unitStates(land[["start"]], transitions), transitions, stageRules = input)

  target <- input("y2000", "opt", land[["initial"]], NULL)[["targetArea"]]
  expect_lt(maxAbs(res[["stateHr"]][, , "urban"], target), 1e-8)
  expect_lt(maxAbs(unitSums(res[["stateHr"]], land[["cluster"]]), unitStates(land[["start"]], transitions)), 1e-8)
  expect_true(all(res[["stateHr"]] >= 0))
  expect_gt(nrow(res[["transitionChange"]]), 0)
  expect_lt(maxAbs(dimSums(res[["stateHr"]], dim = 3), dimSums(land[["initial"]], dim = 3)), 1e-8)
})

# state at the end of one stage "opt" with a target area for urban land (one unit A)
targetStage <- function(state, total, targetArea, allowed) {
  transitions <- oneStage(total)
  rules <- function(year, stage, stateHr, previous, layer) {
    return(list(targetArea = new.magpie(getCells(state), NULL, "urban", fill = targetArea),
                allowedTransitions = allowed))
  }
  cluster <- rep("A", ncells(state))
  res <- toolAllocateTransitions(state, cellMapping(cluster), unitStates(unitSums(state, cluster), transitions),
                                 transitions, rules)
  return(res[["stateHr"]])
}

test_that("target inflow moves back to the giving type where compensation is not possible", {
  # cell 1 has forest and other land, cell 2 only other land; the model takes urban land from
  # forest, the target places it in cell 1, from forest and other land in proportion to their area;
  # forest cannot go to other land, so the target takes forest instead of other land in cell 1
  state <- new.magpie(c("c1", "c2"), "y1995", c("forest", "other", "urban"), fill = c(0.5, 0, 0.5, 1, 0, 0))
  total <- new.magpie("A", NULL, "forest.urban", fill = 0.1)
  res <- targetStage(state, total, c(0.1, 0), c("forest.urban", "other.urban", "other.forest"))
  expect_equal(as.vector(res[, , "urban"]), c(0.1, 0), tolerance = 1e-10)
  expect_equal(as.vector(dimSums(res, dim = 1)), c(0.4, 1.5, 0.1), tolerance = 1e-10)
})

test_that("a cell keeps the target type where no type in the cell can take it", {
  # urban shrinks in cell 1 by the target, but cell 1 holds only urban land and forest, which cannot
  # gain; the model has no outflow of urban land
  state <- new.magpie(c("c1", "c2"), "y1995", c("forest", "other", "urban"), fill = c(0.5, 0, 0, 1, 0.5, 0))
  total <- new.magpie("A", NULL, "other.urban", fill = 0.1)
  allowed <- c("other.urban", "forest.other")
  res <- targetStage(state, total, c(0.3, 0.1), allowed)
  expect_equal(as.vector(dimSums(res, dim = 1)), c(0.5, 0.9, 0.6), tolerance = 1e-10)
  expect_equal(as.vector(res[, , "urban"]), c(0.5, 0.1), tolerance = 1e-10)

  # with a tiny area of a type that can gain, the cell loses the target type to it, and the cluster
  # states still hold
  state["c1", , "other"] <- 1e-33
  res <- expect_silent(targetStage(state, total, c(0.3, 0.1), allowed))
  expect_equal(as.vector(dimSums(res, dim = 1)), c(0.5, 0.9, 0.6), tolerance = 1e-10)
})

test_that("compensation respects allowed", {
  land <- makeLand()
  total <- makeTargetTotals(land[["start"]])
  transitions <- oneStage(total)
  states <- unitStates(land[["start"]], transitions)
  target <- urbanTarget(total, land[["cluster"]])
  allPairs <- paste(rep(landTypes, each = length(landTypes)), rep(landTypes, length(landTypes)), sep = ".")
  run <- function(allowed) {
    input <- function(year, stage, state, previous, layer) {
      out <- target(year, stage, state, previous, layer)
      if (!is.null(out)) {
        out[["allowedTransitions"]] <- allowed
      }
      return(out)
    }
    toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]), states, transitions,
                            stageRules = input)
  }
  area <- function(res, transition) {
    item <- paste0("opt.", transition)
    if (item %in% getNames(res[["transitionsHr"]])) {
      return(sum(res[["transitionsHr"]][, , item]))
    }
    return(0)
  }

  # crop.forest is added by the balancing if it is allowed
  expect_gt(area(run(allPairs), "crop.forest"), 1e-3)
  res <- run(setdiff(allPairs, "crop.forest"))
  expect_equal(area(res, "crop.forest"), 0)
  expect_lt(maxAbs(unitSums(res[["stateHr"]], land[["cluster"]]), states), 1e-8)
  expect_gt(nrow(res[["transitionChange"]]), 0)

  # without crop.forest, crop.other and past.other the other transitions cannot compensate; the
  # target then takes the giving types instead, and the cluster states hold
  forbidden <- c("crop.forest", "crop.other", "past.other")
  expect_silent(res <- run(setdiff(allPairs, forbidden)))
  for (transition in forbidden) {
    expect_equal(area(res, transition), 0)
  }
  expect_lt(maxAbs(unitSums(res[["stateHr"]], land[["cluster"]]), states), 1e-8)
})

test_that("input errors", {
  land <- makeLand()
  good <- oneStage(makeTotals(land[["start"]], c(crop.other = 0.2)))
  states <- unitStates(land[["start"]], good)
  run <- function(initial = land[["initial"]], statesLr = states, transitionsLr = good) {
    toolAllocateTransitions(initial, cellMapping(land[["cluster"]]), statesLr, transitionsLr)
  }
  expect_silent(run())

  bad <- new.magpie(getCells(land[["initial"]]), "y1995", paste0(landTypes, ".a"), fill = 0.1)
  expect_error(run(initial = bad), "must not contain '[.]'")
  bad <- land[["initial"]]
  bad["c1", , "crop"] <- -0.1
  expect_error(run(initial = bad), "startHr must not contain NA or negative values")

  expect_error(run(statesLr = states[, , -1]), "statesLr must be a magclass object with the types of startHr")
  expect_error(run(statesLr = setYears(states, "y2010")), "the units and years of statesLr")
  bad <- good
  getNames(bad) <- "opt.crop.water"
  expect_error(run(transitionsLr = bad), "with types of startHr")
  bad <- new.magpie(c("A", "B", "X"), "y2000", getNames(good), fill = 0.1)
  expect_error(run(transitionsLr = bad), "the units and years of statesLr")

  runWith <- function(constraint) {
    toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]), states, good,
                            constraints = list(a = constraint))
  }
  bound <- land[["initial"]][, , "crop"]
  expect_error(runWith(list(types = "water", upper = bound)), "contains types not in startHr")
  expect_error(runWith(list(types = "crop", upper = bound, lower = bound)), "exactly one of upper and lower")
  expect_error(runWith(list(types = "crop", upper = rep(1, 12))), "bounds must be magclass objects")
  expect_error(runWith(list(types = "crop", upper = land[["initial"]][, , 1:2])), "bounds must be magclass objects")
})

test_that("constraints use the bound of the latest year, and a lower bound limits the other types", {
  startHr <- new.magpie(c("c1", "c2"), "y1995", c("crop", "other"), fill = c(0, 0, 10, 10))
  startLr <- unitSums(startHr, c("A", "A"))
  expansion <- new.magpie("A", "y2015", "opt.other.crop", fill = 2)
  # cropland bound by year: y2015 uses y2010, years before y2000 use y2000
  bound <- new.magpie(c("c1", "c2"), c("y2000", "y2010", "y2020"), "crop", fill = c(10, 10, 0.5, 10, 10, 0.5))
  upper <- list(crop = list(types = "crop", upper = bound))
  res <- toolAllocateTransitions(startHr, cellMapping(c("A", "A")), unitStates(startLr, expansion), expansion,
                                 constraints = upper)
  expect_equal(as.vector(res[["stateHr"]][, , "crop"]), c(0.5, 1.5), tolerance = 1e-6)
  early <- setYears(expansion, "y1995")
  res <- toolAllocateTransitions(startHr, cellMapping(c("A", "A")), unitStates(startLr, early), early,
                                 constraints = upper)
  expect_equal(as.vector(res[["stateHr"]][, , "crop"]), c(1, 1), tolerance = 1e-6)
  # a lower bound on other land (without year) keeps the cropland of cell 1 within the rest of the cell
  lower <- list(other = list(types = "other", lower = new.magpie(c("c1", "c2"), NULL, "other", fill = c(9.5, 0))))
  res <- toolAllocateTransitions(startHr, cellMapping(c("A", "A")), unitStates(startLr, expansion), expansion,
                                 constraints = lower)
  expect_equal(as.vector(res[["stateHr"]][, , "crop"]), c(0.5, 1.5), tolerance = 1e-6)
})

test_that("the deviation from the unit states is reported, and a target is scaled to the unit totals", {
  land <- makeLand()
  transitions <- makeTransitions(land[["start"]])
  states <- unitStates(land[["start"]], transitions)
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]), states, transitions)
  expect_lt(res[["deviationLr"]], 1e-8)
  # unit states that differ from the transitions by 0.01 show as deviation
  shifted <- states
  shifted[, "y2010", "crop"] <- shifted[, "y2010", "crop"] + 0.01
  res <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]), shifted, transitions)
  expect_equal(res[["deviationLr"]], 0.01, tolerance = 1e-6)

  # a target given as distribution (three times the area) gives the same result as the target in Mha
  total <- makeTargetTotals(land[["start"]])
  stages <- oneStage(total)
  rules <- urbanTarget(total, land[["cluster"]])
  scaled <- function(year, stage, state, previous, layer) {
    input <- rules(year, stage, state, previous, layer)
    if (!is.null(input)) {
      input[["targetArea"]] <- 3 * input[["targetArea"]]
    }
    return(input)
  }
  inMha <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                   unitStates(land[["start"]], stages), stages, rules)
  asShare <- toolAllocateTransitions(land[["initial"]], cellMapping(land[["cluster"]]),
                                     unitStates(land[["start"]], stages), stages, scaled)
  expect_lt(maxAbs(inMha[["stateHr"]], asShare[["stateHr"]]), 1e-8)
})

test_that("bounds are matched by grid cell", {
  startHr <- new.magpie(c("c1", "c2"), "y1995", c("crop", "other"), fill = c(0, 0, 10, 10))
  expansion <- new.magpie("A", "y2015", "opt.other.crop", fill = 2)
  states <- unitStates(unitSums(startHr, c("A", "A")), expansion)
  .run <- function(bound) {
    toolAllocateTransitions(startHr, cellMapping(c("A", "A")), states, expansion,
                            constraints = list(crop = list(types = "crop", upper = bound)))[["stateHr"]]
  }
  ordered <- .run(new.magpie(c("c1", "c2"), "y2010", "crop", fill = c(0.5, 10)))
  reversed <- .run(new.magpie(c("c2", "c1"), "y2010", "crop", fill = c(10, 0.5)))
  expect_equal(reversed, ordered)
  expect_equal(as.vector(ordered[, , "crop"]), c(0.5, 1.5), tolerance = 1e-6)
})
