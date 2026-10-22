#' calcLandTransitionsInput
#'
#' Land transitions of an input dataset for the downscaling along its transitions
#' (\code{\link{calcLandTransitionsHighRes}}): the transitions of each time step by stage, in the order
#' in which the model applies them, the unit states, the start map and the rules that place the stages
#' on the grid cells.
#'
#' @details
#' input = "magpie": the cluster transitions by stage and sub-type of \code{magpie4::landTransitions}
#' (forestry and other land split into the sub-types of modules 32 and 35). The start map is the land-use
#' initialisation of the run (\code{avl_land_full_t_0.5.mz}), with forestry and other land split into
#' equal parts of their sub-types (\code{\link{toolAllocateTransitions}} fits them to the cluster start).
#'
#' Rules (see \code{\link{calcLandTransitionsHighRes}}):
#' \itemize{
#'   \item recovery is placed on the othernat inflow of the previous year (snvRelocation, peatRewetting
#'   and optimisation), recoveryForestry on the inflow from forestry;
#'   \item the SNV relocation in proportion to the relocation target (\code{SNVTargetCropland_0.5.mz},
#'   interpolated for the SNV share of the run as \code{i29_snv_relocation_target}) not relocated in
#'   earlier years, if the run has the target in its output folder;
#'   \item the land reduction on drained peat in proportion to the drained peat of the grid cells
#'   (peatland, updated each year by \code{\link{toolDownscalePeatland}});
#'   \item in the optimisation, urban land follows the 0.5 degree data of the run
#'   (\code{c34_urban_scenario}); the other transitions are changed within the rules of the land module
#'   and without new inflow to secondary forest (only harvest regrowth and restoration), youngsecdf and
#'   primary forest.
#' }
#'
#' @param input name of an input dataset, currently only "magpie"
#' @return List with
#' \item{mainOf}{named character vector: the land type of each sub-type (the type itself for land types
#' without sub-types)}
#' \item{mapping}{data frame with columns cell and cluster}
#' \item{startHr}{magclass object of the grid cell state before the first year (grid cells x one year x
#' sub-types)}
#' \item{statesLr}{magclass object of the unit states at the end of each year (units x years x sub-types)}
#' \item{transitionsLr}{magclass object of the unit transitions (units x years x "stage.from.to")}
#' \item{splitLr}{named list of the unit compositions of items without transitions (magclass objects,
#' units x years x parts), applied in order: cropland into crop area, fallow and tree cover, then crop
#' area into crops by water supply}
#' \item{rules}{list of rules, each a list with stage, kind and its data}
#' @author Patrick v. Jeetze
#' @seealso \code{\link{calcLandConstraints}}, \code{\link{calcLandTransitionsHighRes}}
calcLandTransitionsInput <- function(input) {
  # Initial grid cell state by sub-type: pasture and rangeland merged, forestry and other land split into
  # equal parts; toolAllocateTransitions then fits each sub-type to its cluster start total.
  .initialSubtypes <- function(initial, mainOf) {
    forestryTypes <- names(mainOf)[mainOf == "forestry"]
    otherTypes <- names(mainOf)[mainOf == "other"]
    initial[initial < 0] <- 0
    out <- new.magpie(getCells(initial), getYears(initial), names(mainOf), fill = 0)
    for (k in c("crop", "primforest", "secdforest", "urban")) {
      out[, , k] <- initial[, , k]
    }
    out[, , "past"] <- dimSums(initial[, , c("past", "range")], dim = 3)
    for (k in forestryTypes) {
      out[, , k] <- initial[, , "forestry"] / length(forestryTypes)
    }
    for (k in otherTypes) {
      out[, , k] <- dimSums(initial[, , c("primother", "secdother")], dim = 3) / length(otherTypes)
    }
    return(out)
  }

  # The same items with the data dimension set "type", for binding objects of different sets.
  .asTypes <- function(x) {
    getSets(x, fulldim = FALSE)[3] <- "type"
    return(x)
  }

  if (input == "magpie") {
    # =================================================================================================
    # Inputs of the run
    # =================================================================================================

    cfg <- readSource("MagpieRunFiles", subtype = "config")
    mapping <- readSource("MagpieFulldataGdx", subtype = "clustermap")
    dataLr <- readSource("MagpieFulldataGdx", subtype = "landAccounting")
    transitionsLr <- readSource("MagpieFulldataGdx", subtype = "landTransitions")
    landStartHr <- readSource("MagpieRunFiles", subtype = "landuseInit")[, "y1995", ]

    clusters <- getItems(transitionsLr, dim = 1)
    cells <- getCells(landStartHr)
    cluster <- as.character(mapping[["cluster"]][match(cells, mapping[["cell"]])])
    stopifnot("mapping lacks grid cells of landStartHr" = !anyNA(cluster),
              "clusters of the transitions and the mapping differ" = setequal(clusters, unique(cluster)))
    years <- getYears(transitionsLr)

    # =================================================================================================
    # Land types and cluster states
    # =================================================================================================

    # main land type of each sub-type of the run: forestry types of module 32 and other land types of
    # module 35; the other land types are their own sub-type
    forestryTypes <- getItems(dataLr[["forestry"]], "type32")
    otherTypes <- getItems(dataLr[["other"]], "othertype35")
    mainOf <- c(crop = "crop", past = "past",
                stats::setNames(rep("forestry", length(forestryTypes)), forestryTypes),
                primforest = "primforest", secdforest = "secdforest", urban = "urban",
                stats::setNames(rep("other", length(otherTypes)), otherTypes))
    types <- names(mainOf)

    # cluster state at the end of each year by sub-type
    statesLr <- mbind(.asTypes(dataLr[["land"]][, years, c("crop", "past", "primforest", "secdforest", "urban")]),
                      .asTypes(dataLr[["forestry"]][, years, ]),
                      .asTypes(dimSums(dataLr[["other"]][, years, ], dim = "ac")))[, , types]

    # composition of cropland in each cluster: crop area, fallow and tree cover, and crop area by crop
    # and water supply
    cropParts <- c("crop_area", "crop_fallow", "crop_treecover")
    splitLr <- list(crop = readSource("MagpieFulldataGdx", subtype = "land")[, , cropParts],
                    crop_area = readSource("MagpieFulldataGdx", subtype = "crop"))

    # =================================================================================================
    # Rules
    # =================================================================================================

    # recovery on the othernat inflow of the previous year, forestry abandonment first
    inflowStages <- c("snvRelocation", "peatRewetting", "optimisation")
    rules <- list(
      list(stage = "recoveryForestry", kind = "previousInflow", type = "othernat", fromStages = inflowStages,
           fromTypes = forestryTypes),
      list(stage = "recovery", kind = "previousInflow", type = "othernat", fromStages = inflowStages)
    )

    # SNV relocation in proportion to the cropland that requires relocation under the SNV policy
    # (Copernicus land cover), interpolated between the 20 % and 50 % targets for the SNV share of the
    # grid cell as i29_snv_relocation_target; runs without these data at 0.5 degree in their output folder
    # allocate the relocation like the other transitions, within the crop constraint
    stageOf <- getItems(transitionsLr, dim = 3.1, full = TRUE)
    snvRelocation <- "snvRelocation" %in% stageOf && sum(transitionsLr[, , "snvRelocation"]) > 0
    if (snvRelocation && !"snvTarget" %in% readSource("MagpieRunFiles", subtype = "available")) {
      warning("The run has no SNV target cropland at 0.5 degree (SNVTargetCropland_0.5.mz): ",
              "the SNV relocation is allocated without it")
    } else if (snvRelocation) {
      snv <- readSource("MagpieFulldataGdx", subtype = "snv")
      targetData <- readSource("MagpieRunFiles", subtype = "snvTarget")
      stopifnot("SNVTargetCropland_0.5.mz lacks grid cells of the start map" = all(cells %in% getCells(targetData)))
      targetData <- setYears(targetData[cells, , ], NULL)
      share <- new.magpie(cells, NULL, fill = unname(snv[["shareIso"]][getItems(targetData, dim = 1.3, full = TRUE)]))
      target20 <- collapseDim(targetData[, , "SNV20TargetCropland"], dim = 3)
      target50 <- collapseDim(targetData[, , "SNV50TargetCropland"], dim = 3)
      x1 <- snv[["relocationX1"]]
      x2 <- snv[["relocationX2"]]
      target <- pmin(share / x1, 1) * target20 + pmax(share - x1, 0) / (x2 - x1) * (target50 - target20)
      rules[[length(rules) + 1]] <- list(stage = "snvRelocation", kind = "depletingTarget", type = "crop",
                                         targetHr = target)
    }

    # land reduction on drained peat in proportion to the drained peat of the grid cells (peatland realisation
    # v2), from the initial grid cell peat by category; categories without data start at 0
    if (cfg[["gms"]][["peatland"]] == "v2") {
      peatlandLr <- readSource("MagpieFulldataGdx", subtype = "peatland")
      peatlandStartHr <- readSource("MagpieRunFiles", subtype = "peatland")
      stopifnot("peatland data lack grid cells of the start map" = all(cells %in% getCells(peatlandStartHr)))
      categories <- getNames(peatlandLr)
      peatStartHr <- new.magpie(cells, NULL, categories, fill = 0)
      for (k in intersect(categories, getNames(peatlandStartHr))) {
        peatStartHr[, , k] <- setYears(peatlandStartHr[cells, 1, k], NULL)
      }
      rules[[length(rules) + 1]] <- list(stage = "peatRewetting", kind = "peatland", startHr = peatStartHr,
                                         lowRes = peatlandLr,
                                         onTypes = intersect(categories, c("crop", "past", "forestry")))
    }

    # urban land from the data (scaled to the cluster urban land of the model); compensation within the
    # transitions with a positive upper bound in the land module, and no new inflow to secondary forest
    # (only harvest regrowth and restoration), youngsecdf and primary forest
    if (cfg[["gms"]][["urban"]] == "exo_nov21") {
      urbanHr <- readSource("MagpieRunFiles", subtype = "urban")
      stopifnot("urban data lack grid cells of the start map" = all(cells %in% getCells(urbanHr)))
      urbanHr <- urbanHr[cells, , cfg[["gms"]][["c34_urban_scenario"]]]
      getNames(urbanHr) <- "urban"
      upper <- dataLr[["transitionsUpper"]]
      mainAllowed <- getNames(upper)[as.vector(dimSums((upper > 0) * 1, dim = c(1, 2)) > 0)]
      pairs <- expand.grid(from = types, to = types, stringsAsFactors = FALSE)
      allowed <- pairs[["from"]] != pairs[["to"]] &
        paste(mainOf[pairs[["from"]]], mainOf[pairs[["to"]]], sep = ".") %in% mainAllowed &
        !pairs[["to"]] %in% c("secdforest", "youngsecdf", "primforest")
      allowedTransitions <- paste(pairs[["from"]], pairs[["to"]], sep = ".")[allowed]
      rules[[length(rules) + 1]] <- list(stage = "optimisation", kind = "targetArea", targetHr = urbanHr,
                                         allowedTransitions = allowedTransitions)
    }

    out <- list(mainOf = mainOf, mapping = mapping, startHr = .initialSubtypes(landStartHr, mainOf),
                statesLr = statesLr, transitionsLr = transitionsLr, splitLr = splitLr, rules = rules)
  } else {
    stop("Unsupported input dataset \"", input, "\"")
  }

  return(list(x = out, class = "list", unit = "Mha",
              description = "Land transitions of the input by year and stage, unit states, start map and rules"))
}
