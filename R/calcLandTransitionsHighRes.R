#' calcLandTransitionsHighRes
#'
#' Land of an input dataset on the grid along its transitions: the unit transitions of each time step
#' (\code{\link{calcLandTransitionsHarmonized}}) are allocated to the grid cells stage by stage with
#' \code{\link{toolAllocateTransitions}}, within the constraints of \code{\link{calcLandConstraints}}.
#' The grid cell transitions sum exactly to the unit transitions, so the grid cell states sum to the
#' unit states.
#'
#' @details
#' The input declares the rules that place a stage, at most one per stage:
#' \itemize{
#'   \item previousInflow: the outflow of \code{type} in proportion to its inflow in the previous year
#'   (from the stages \code{fromStages}, and only from \code{fromTypes} if given) per area of the type;
#'   \item depletingTarget: the outflow of \code{type} in proportion to the part of the grid cell target
#'   \code{targetHr} not yet taken by this stage in earlier years, per area of the type;
#'   \item peatland: the outflow of the land types with drained peat (\code{onTypes}) in proportion to
#'   the drained peat of the previous year per area of the land type; the peat (start \code{startHr},
#'   unit values \code{lowRes}) follows the land year by year with \code{\link{toolDownscalePeatland}}, its
#'   decline on the land types taken from the stage of the rule;
#'   \item targetArea: the grid cell distribution of a type at the end of the stage (\code{targetHr},
#'   grid cells x years), with the transitions that the compensation may add
#'   (\code{allowedTransitions}), see \code{\link{toolAllocateTransitions}}.
#' }
#' Grid cells without inflow, target or peat take the rest (weight 1e-6) in proportion to the area of
#' the type. Items without transitions
#' (\code{splitLr}) are split with the unit shares after the allocation: every grid cell of a unit gets
#' the composition of the unit; units without area take the composition of all units, and the start
#' takes the composition of the first year.
#'
#' @inheritParams calcLandTransitionsHarmonized
#' @param harmonizationPeriod Two integer values, see \code{\link{calcLandHarmonized}}; not used with
#' \code{harmonization = "none"}, where callers leave it out so that they share one cached result
#' @return List with
#' \item{landHr}{magclass object of grid cell land (Mha) by land type and sub-type (names "land.type", type
#' "total" for land types without sub-types), the state before the first year labelled y1985}
#' \item{splitHr}{named list of the items without transitions split by the unit shares (magclass objects,
#' Mha), as \code{splitLr} of \code{\link{calcLandTransitionsInput}}}
#' \item{transitionsHr}{magclass object of the gross grid cell transitions per year (Mha) by stage, land
#' type and sub-type (names "stage.from.fromType.to.toType", stages in the order of the input)}
#' \item{diagnostics, excess, transitionChange, deviationLr}{as returned by
#' \code{\link{toolAllocateTransitions}}}
#' \item{peatlandHr, peatlandExcess}{grid cell peatland by category (Mha) and drained peat above the grid
#' cell land, or \code{NULL} without a peatland rule}
#' @author Patrick v. Jeetze
#' @seealso \code{\link{calcLandHighRes}}, \code{\link{fullDOWNSCALEDMAGPIE}}
calcLandTransitionsHighRes <- function(input, target = "landuseinit", harmonizationPeriod = NULL, # nolint: cyclocomp_linter
                                       harmonization = "none") {
  # =================================================================================================
  # Rules
  # =================================================================================================

  # Weight of an outflow in proportion to amount per area; grid cells without amount take the rest (1e-6)
  # in proportion to their area.
  .perArea <- function(amount, area) {
    weight <- amount / area
    weight[area <= 0] <- 0
    return(weight + 1e-6)
  }

  # Weight of the outflow of rule[["type"]] by its inflow in the previous year per area.
  .previousInflowWeight <- function(rule, state, previous) {
    select <- getItems(previous, dim = 3.1, full = TRUE) %in% rule[["fromStages"]] &
      getItems(previous, dim = 3.3, full = TRUE) == rule[["type"]]
    if (!is.null(rule[["fromTypes"]])) {
      select <- select & getItems(previous, dim = 3.2, full = TRUE) %in% rule[["fromTypes"]]
    }
    inflow <- 0
    if (any(select)) {
      inflow <- setYears(dimSums(previous[, , select], dim = 3), NULL)
    }
    return(.perArea(inflow, state[, , rule[["type"]]]))
  }

  # Grid cell area of land types, summed over their sub-types.
  .mainLand <- function(state, onTypes) {
    land <- list()
    for (k in onTypes) {
      land[[k]] <- setNames(dimSums(state[, , names(mainOf)[mainOf == k]], dim = 3), k)
    }
    return(mbind(land))
  }

  # Weight of the outflow of the sub-types of the land types with drained peat: the drained peat per area
  # of the land type in the grid cell.
  .peatlandWeight <- function(state, peat, onTypes) {
    weight <- .perArea(peat[, , onTypes], .mainLand(state, onTypes))
    subTypes <- names(mainOf)[mainOf %in% onTypes]
    out <- state[, , subTypes] * 0
    for (k in subTypes) {
      out[, , k] <- collapseDim(weight[, , mainOf[[k]]], dim = 3)
    }
    return(out)
  }

  # Labels "land.type" of the sub-types (type "total" for land types without sub-types).
  .typeLabels <- function() {
    subType <- ifelse(names(mainOf) == mainOf, "total", names(mainOf))
    return(stats::setNames(paste(mainOf, subType, sep = "."), names(mainOf)))
  }

  # Shares of the items in each unit and year; units without area take the shares of all units, copied
  # to the grid cells of each unit, the first year also used for the start (y1985).
  .gridShares <- function(x, mapping) {
    filled <- x + (dimSums(x, dim = 3) == 0) * dimSums(x, dim = 1)
    shares <- filled / dimSums(filled, dim = 3)
    shares[is.nan(shares)] <- 0
    return(toolAggregate(mbind(setYears(shares[, 1, ], "y1985"), shares), mapping, from = "cluster", to = "cell"))
  }

  # =================================================================================================
  # Inputs
  # =================================================================================================

  transitions <- calcOutput("LandTransitionsHarmonized", input = input, target = target,
                            harmonizationPeriod = harmonizationPeriod, harmonization = harmonization,
                            aggregate = FALSE)
  constraints <- calcOutput("LandConstraints", input = input, aggregate = FALSE)
  mainOf <- transitions[["mainOf"]]
  mapping <- transitions[["mapping"]]

  # main types in the constraints stand for their sub-types
  constraints <- lapply(constraints, function(constraint) {
    types <- constraint[["types"]]
    constraint[["types"]] <- unique(c(setdiff(types, mainOf), names(mainOf)[mainOf %in% types]))
    return(constraint)
  })

  # rules by stage, and the layer of toolAllocateTransitions: the area each depleting target stage has
  # taken so far, and the peat of the grid cells
  ruleOf <- list()
  layer <- list(used = list(), peatland = NULL)
  for (rule in transitions[["rules"]]) {
    stopifnot("at most one rule per stage" = is.null(ruleOf[[rule[["stage"]]]]))
    ruleOf[[rule[["stage"]]]] <- rule
    if (rule[["kind"]] == "depletingTarget") {
      layer[["used"]][[rule[["stage"]]]] <- 0
    } else if (rule[["kind"]] == "peatland") {
      peatlandRule <- rule
      layer[["peatland"]] <- list(current = rule[["startHr"]], out = list(), excess = list())
    }
  }

  # =================================================================================================
  # Allocation
  # =================================================================================================

  stageRules <- function(year, stage, state, previous, layer) {
    stageInput <- list()
    rule <- ruleOf[[stage]]
    if (is.null(rule)) {
      return(stageInput)
    }
    if (rule[["kind"]] == "previousInflow") {
      if (!is.null(previous)) {
        stageInput[["fromWeight"]] <- .previousInflowWeight(rule, state, previous)
      }
    } else if (rule[["kind"]] == "depletingTarget") {
      type <- rule[["type"]]
      stageInput[["fromWeight"]] <- .perArea(pmax(rule[["targetHr"]] - layer[["used"]][[stage]], 0), state[, , type])
    } else if (rule[["kind"]] == "peatland") {
      stageInput[["fromWeight"]] <- .peatlandWeight(state, layer[["peatland"]][["current"]], rule[["onTypes"]])
    } else if (rule[["kind"]] == "targetArea") {
      stageInput[["targetArea"]] <- rule[["targetHr"]][, year, ]
      stageInput[["allowedTransitions"]] <- rule[["allowedTransitions"]]
    } else {
      stop("Unsupported rule kind \"", rule[["kind"]], "\"")
    }
    return(stageInput)
  }

  # Updates the layer at the end of each year: the area taken by the depleting target stages, and the peat
  # on the land and the land reduction of the stage of the peatland rule.
  updateLayer <- function(year, state, flows, layer) {
    stageOf <- character(0)
    fromMain <- character(0)
    if (!is.null(flows)) {
      stageOf <- getItems(flows, dim = 3.1, full = TRUE)
      fromMain <- mainOf[getItems(flows, dim = 3.2, full = TRUE)]
    }
    for (stage in names(layer[["used"]])) {
      if (stage %in% stageOf) {
        layer[["used"]][[stage]] <- layer[["used"]][[stage]] +
          setYears(dimSums(flows[, , stageOf == stage], dim = 3), NULL)
      }
    }

    if (!is.null(layer[["peatland"]])) {
      onTypes <- peatlandRule[["onTypes"]]
      landHr <- .mainLand(state, onTypes)
      # land reduction on drained peat by land type
      loss <- landHr * 0
      for (k in onTypes) {
        reduction <- stageOf == peatlandRule[["stage"]] & fromMain == k
        if (any(reduction)) {
          loss[, , k] <- dimSums(flows[, , reduction], dim = 3)
        }
      }

      peat <- layer[["peatland"]]
      result <- toolDownscalePeatland(landHr, peatlandRule[["lowRes"]][, year, ], setYears(peat[["current"]], year),
                                      mapping, loss)
      peat[["current"]] <- setYears(result[["peatlandHr"]], NULL)
      peat[["out"]][[year]] <- result[["peatlandHr"]]
      peat[["excess"]][[year]] <- result[["excess"]]
      layer[["peatland"]] <- peat
    }
    return(layer)
  }

  # the unit totals hold exactly (see .fitSourceArea in toolAllocateTransitions), so the constraints
  # need not converge further than 1e-6 Mha per unit
  result <- toolAllocateTransitions(transitions[["startHr"]], mapping, transitions[["statesLr"]],
                                    transitions[["transitionsLr"]], stageRules, constraints, layer = layer,
                                    updateLayer = updateLayer, tol = 1e-6)

  # =================================================================================================
  # Outputs
  # =================================================================================================

  # grid cell land by land type and sub-type ("land.type"), the start labelled y1985
  label <- .typeLabels()
  landHr <- mbind(setYears(result[["startHr"]], "y1985"), result[["stateHr"]])
  getNames(landHr) <- unname(label[getNames(landHr)])
  getSets(landHr, fulldim = FALSE) <- c("x.y.iso", "year", "land.type")

  # gross grid cell transitions by stage, land type and sub-type
  transitionsHr <- result[["transitionsHr"]]
  # the largest output: without the second reference, the renaming does not keep two copies of it
  result[["transitionsHr"]] <- NULL
  getNames(transitionsHr) <- paste(getItems(transitionsHr, dim = 3.1, full = TRUE),
                                   label[getItems(transitionsHr, dim = 3.2, full = TRUE)],
                                   label[getItems(transitionsHr, dim = 3.3, full = TRUE)], sep = ".")
  getSets(transitionsHr, fulldim = FALSE) <- c("x.y.iso", "year", "stage.from.fromType.to.toType")

  # items without transitions split with the unit shares; the first split divides a land type, later
  # ones an item of an earlier split
  splitHr <- list()
  for (item in names(transitions[["splitLr"]])) {
    if (item %in% mainOf) {
      parent <- setNames(dimSums(landHr[, , item], dim = 3), NULL)
    } else {
      parent <- NULL
      for (earlier in names(splitHr)) {
        if (item %in% getItems(splitHr[[earlier]], dim = 3)) {
          parent <- setNames(splitHr[[earlier]][, , item], NULL)
        }
      }
      stopifnot("a split divides a land type or an item of an earlier split" = !is.null(parent))
    }
    splitHr[[item]] <- .gridShares(transitions[["splitLr"]][[item]], mapping) * parent
  }

  peatlandHr <- NULL
  peatlandExcess <- NULL
  peat <- result[["layer"]][["peatland"]]
  if (!is.null(peat)) {
    peatlandHr <- mbind(peat[["out"]])
    peatlandExcess <- do.call(rbind, unname(peat[["excess"]]))
  }

  out <- list(landHr = landHr, splitHr = splitHr, transitionsHr = transitionsHr,
              diagnostics = result[["diagnostics"]], excess = result[["excess"]],
              transitionChange = result[["transitionChange"]], deviationLr = result[["deviationLr"]],
              peatlandHr = peatlandHr, peatlandExcess = peatlandExcess)
  toolExpectTrue(out[["deviationLr"]] <= 1e-5,
                 paste0("Unit states equal the states of the input (largest deviation ",
                        signif(out[["deviationLr"]], 3), " Mha)"))

  return(list(x = out, class = "list", unit = "Mha",
              description = "Land, transitions and peatland of the input downscaled along its transitions"))
}
