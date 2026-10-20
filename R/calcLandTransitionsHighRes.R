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
    return(ifelse(area > 0, amount / area, 0) + 1e-6)
  }

  # Weight of the outflow of rule[["type"]] by its inflow in the previous year per area.
  .previousInflowWeight <- function(rule, state, previousFlows) {
    inflow <- 0
    for (flows in previousFlows[rule[["fromStages"]]]) {
      if (is.null(flows)) {
        next
      }
      select <- sub("^.*[.]", "", colnames(flows)) == rule[["type"]]
      if (!is.null(rule[["fromTypes"]])) {
        select <- select & sub("[.].*$", "", colnames(flows)) %in% rule[["fromTypes"]]
      }
      inflow <- inflow + rowSums(flows[, select, drop = FALSE])
    }
    weight <- .perArea(inflow, state[, rule[["type"]]])
    return(matrix(weight, ncol = 1, dimnames = list(NULL, rule[["type"]])))
  }

  # Grid cell area of land types, summed over their sub-types (grid cells x types).
  .mainLand <- function(state, onTypes, mainOf) {
    land <- matrix(0, nrow(state), length(onTypes), dimnames = list(rownames(state), onTypes))
    for (k in onTypes) {
      land[, k] <- rowSums(state[, names(mainOf)[mainOf == k], drop = FALSE])
    }
    return(land)
  }

  # Weight of the outflow of the sub-types of the land types with drained peat: the drained peat per area
  # of the land type in the grid cell (peat: grid cells x categories).
  .peatlandWeight <- function(state, peat, onTypes, mainOf) {
    weight <- .perArea(peat[, onTypes, drop = FALSE], .mainLand(state, onTypes, mainOf))
    subTypes <- names(mainOf)[mainOf %in% onTypes]
    weight <- weight[, mainOf[subTypes], drop = FALSE]
    colnames(weight) <- subTypes
    return(weight)
  }

  # =================================================================================================
  # Outputs
  # =================================================================================================

  # Grid cell state by sub-type as magclass object with names "land.type" (type "total" for land types
  # without sub-types), the start labelled y1985.
  .landOutput <- function(state, mainOf, cells) {
    out <- aperm(state, c(1, 3, 2))
    dimnames(out) <- list(cells, c("y1985", dimnames(state)[[3]][-1]), unname(.typeLabels(mainOf)))
    return(.asCellMagpie(out, "land.type"))
  }

  # Gross grid cell transitions per year by stage and sub-type, as magclass object with the names
  # "stage.from.fromType.to.toType" in the order of the stages; transitions without flow in any grid cell
  # and year are left out.
  .transitionOutput <- function(transitions, mainOf, cells) {
    years <- names(transitions)
    label <- .typeLabels(mainOf)

    # flows of each year (grid cells x "stage.from.fromType.to.toType"), without columns of zeros
    byYear <- list()
    for (t in years) {
      byYear[[t]] <- matrix(0, length(cells), 0)
      for (stage in names(transitions[[t]])) {
        flows <- transitions[[t]][[stage]]
        flows <- flows[, colSums(flows) > 0, drop = FALSE]
        if (ncol(flows) == 0) {
          next
        }
        colnames(flows) <- paste(stage, label[sub("[.].*$", "", colnames(flows))],
                                 label[sub("^.*[.]", "", colnames(flows))], sep = ".")
        byYear[[t]] <- cbind(byYear[[t]], flows)
      }
    }

    # items in the order of the stages, so that readers can replay them
    itemNames <- unique(unlist(lapply(byYear, colnames)))
    stageNames <- unique(unlist(lapply(transitions, names)))
    itemNames <- itemNames[order(match(sub("[.].*$", "", itemNames), stageNames))]

    out <- array(0, dim = c(length(cells), length(years), length(itemNames)),
                 dimnames = list(cells, years, itemNames))
    for (t in years) {
      out[, t, colnames(byYear[[t]])] <- byYear[[t]]
    }

    return(.asCellMagpie(out, "stage.from.fromType.to.toType"))
  }

  # Labels "land.type" of the sub-types (type "total" for land types without sub-types).
  .typeLabels <- function(mainOf) {
    subType <- ifelse(names(mainOf) == mainOf, "total", names(mainOf))
    return(stats::setNames(paste(mainOf, subType, sep = "."), names(mainOf)))
  }

  # magclass object of an array grid cells x years x items, with the sets of the 0.5 degree outputs.
  .asCellMagpie <- function(x, itemSet) {
    out <- as.magpie(x, spatial = 1, temporal = 2)
    getSets(out, fulldim = FALSE) <- c("x.y.iso", "year", itemSet)
    return(out)
  }

  # magclass object of one year of a matrix (grid cells x items).
  .asYear <- function(values, year) {
    x <- array(values, c(nrow(values), 1, ncol(values)), list(rownames(values), year, colnames(values)))
    return(.asCellMagpie(x, "data"))
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
  cells <- rownames(transitions[["startHr"]])

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
      weight <- .perArea(pmax(rule[["targetHr"]] - layer[["used"]][[stage]], 0), state[, type])
      stageInput[["fromWeight"]] <- matrix(weight, ncol = 1, dimnames = list(NULL, type))
    } else if (rule[["kind"]] == "peatland") {
      stageInput[["fromWeight"]] <- .peatlandWeight(state, layer[["peatland"]][["current"]], rule[["onTypes"]],
                                                    mainOf)
    } else if (rule[["kind"]] == "targetArea") {
      stageInput[["targetArea"]] <- list()
      for (type in names(rule[["targetHr"]])) {
        stageInput[["targetArea"]][[type]] <- rule[["targetHr"]][[type]][, year]
      }
      stageInput[["allowedTransitions"]] <- rule[["allowedTransitions"]]
    } else {
      stop("Unsupported rule kind \"", rule[["kind"]], "\"")
    }
    return(stageInput)
  }

  # Updates the layer at the end of each year: the area taken by the depleting target stages, and the peat
  # on the land and the land reduction of the stage of the peatland rule.
  updateLayer <- function(year, state, flows, layer) {
    for (stage in names(layer[["used"]])) {
      if (!is.null(flows[[stage]])) {
        layer[["used"]][[stage]] <- layer[["used"]][[stage]] + rowSums(flows[[stage]])
      }
    }

    if (!is.null(layer[["peatland"]])) {
      onTypes <- peatlandRule[["onTypes"]]
      reduction <- flows[[peatlandRule[["stage"]]]]
      fromMain <- mainOf[sub("[.].*$", "", colnames(reduction))]
      loss <- matrix(0, length(cells), length(onTypes), dimnames = list(cells, onTypes))
      for (k in onTypes) {
        loss[, k] <- rowSums(reduction[, fromMain == k, drop = FALSE])
      }

      peat <- layer[["peatland"]]
      result <- toolDownscalePeatland(.asYear(.mainLand(state, onTypes, mainOf), year),
                                      peatlandRule[["lowRes"]][, year, ], .asYear(peat[["current"]], year),
                                      mapping, .asYear(loss, year))
      peat[["current"]][] <- as.array(result[["peatlandHr"]])[, 1, ]
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

  landHr <- .landOutput(result[["stateHr"]], mainOf, cells)

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

  out <- list(landHr = landHr, splitHr = splitHr,
              transitionsHr = .transitionOutput(result[["transitionsHr"]], mainOf, cells),
              diagnostics = result[["diagnostics"]], excess = result[["excess"]],
              transitionChange = result[["transitionChange"]], deviationLr = result[["deviationLr"]],
              peatlandHr = peatlandHr, peatlandExcess = peatlandExcess)
  toolExpectTrue(out[["deviationLr"]] <= 1e-5,
                 paste0("Unit states equal the states of the input (largest deviation ",
                        signif(out[["deviationLr"]], 3), " Mha)"))

  return(list(x = out, class = "list", unit = "Mha",
              description = "Land, transitions and peatland of the input downscaled along its transitions"))
}
