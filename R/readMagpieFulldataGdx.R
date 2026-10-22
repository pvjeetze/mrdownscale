#' readMagpieFulldataGdx
#'
#' Read function for data coming from the MAgPIE model.
#'
#' @param subtype type of data to be read in. Available options are
#' land, crop, woodHarvestWeight, woodHarvestArea, fertilizerRegional, clustermap, and for the
#' downscaling along the transitions of the model (lists of gdx symbols at cluster level):
#' landAccounting (cluster land by type and sub-type, upper bounds of the land transitions), landTransitions,
#' ageClasses and ageClassTransitions (\code{magpie4::landTransitions}, \code{magpie4::ageClasses},
#' \code{magpie4::ageClassTransitions}), peatland, landConservation, snv, bii, climateFixYear
#' @author Pascal Sauer, Jan Philipp Dietrich, Patrick v. Jeetze
readMagpieFulldataGdx <- function(subtype) {
  "!# @monitor magpie4:::land"
  "!# @monitor magpie4:::croparea"
  "!# @monitor magpie4:::woodProduction"
  "!# @monitor magpie4:::woodHarvestArea"
  "!# @monitor magpie4:::NitrogenBudget"
  "!# @monitor magpie4:::PeatlandArea"
  "!# @monitor magpie4:::BII"
  "!# @monitor magpie4:::ageClasses"
  "!# @monitor magpie4:::landTransitions"
  "!# @monitor magpie4:::landConservation"
  "!# @monitor magpie4:::ageClassTransitions"

  gdx <- "fulldata.gdx"
  stopifnot(file.exists(gdx),
            length(Sys.glob("clustermap_*.rds")) == 1)

  clustermap <- readRDS(Sys.glob("clustermap_*.rds"))

  # Level of a variable of the run.
  .lvl <- function(symbol) {
    return(gdx2::readGDX(gdx, symbol, select = list(type = "level")))
  }

  if (subtype == "land") {
    x <- magpie4::land(gdx, level = "cell", subcategories = c("crop", "forestry"))
    x <- add_columns(x, c("forestry_aff_pltns", "forestry_aff_secdf"), fill = 0)
    if (as.vector(gdx2::readGDX(gdx, "s32_aff_plantation")) == 1) {
      x[, , "forestry_aff_pltns"] <- x[, , "forestry_aff"]
    } else {
      x[, , "forestry_aff_secdf"] <- x[, , "forestry_aff"]
    }
    x <- x[, , "forestry_aff", invert = TRUE]
    getSets(x) <- c("region", "id", "year", "data")
    unit <- "Mha"
    description <- "land use information"
  } else if (subtype == "crop") {
    x <- magpie4::croparea(gdx, level = "cell", product_aggr = FALSE, water_aggr = FALSE)
    attr(x, "gdxdata") <- NULL
    getSets(x) <- c("region", "id", "year", "data", "water")
    unit <- "Mha"
    description <- "crop land information separated by irrigated/rainfed"
  } else if (subtype == "woodHarvestWeight") {
    x <- magpie4::woodProduction(gdx)
    stopifnot(identical(getComment(x), " unit: Pg DM yr-1"))
    unit <- "Pg DM yr-1"
    description <- "roundwood and fuelwood harvest weight separated by source"
  } else if (subtype == "woodHarvestArea") {
    x <- magpie4::woodHarvestArea(gdx)
    stopifnot(identical(getComment(x), " unit: Mha yr-1"))
    unit <- "Mha yr-1"
    description <- "wood harvest area separated by source and age classes"
  } else if (subtype == "fertilizer") {
    x <- magpie4::NitrogenBudget(gdx, level = "cell", cropTypes = TRUE,
                                 threshold = 0.0023, progress = FALSE)
    x <- collapseDim(x[, , "fertilizer"])
    getSets(x) <- c("region", "id", "year", "data")
    unit <- "Tg yr-1"
    description <- "fertilizer per croptype"
  } else if (subtype == "clustermap") {
    return(list(x = clustermap, class = "data.frame"))
  } else if (subtype == "landAccounting") {
    # cluster land by type and by the sub-types of forestry and other land, and the upper bounds of the land
    # transitions
    x <- list(land = .lvl("ov_land"), forestry = .lvl("ov_land_forestry"), other = .lvl("ov_land_other"),
              transitionsUpper = gdx2::readGDX(gdx, "ov_lu_transitions", select = list(type = "upper")))
    return(list(x = x, class = "list"))
  } else if (subtype == "landTransitions") {
    x <- magpie4::landTransitions(gdx, level = "cell")
    unit <- "Mha"
    description <- "land transitions by stage and sub-type in each time step"
  } else if (subtype == "ageClasses") {
    x <- magpie4::ageClasses(gdx, level = "cell")
    unit <- "Mha"
    description <- "area by land type and age class"
  } else if (subtype == "ageClassTransitions") {
    x <- magpie4::ageClassTransitions(gdx, level = "cell")
    unit <- "Mha"
    description <- "flows by process, land type and age class in each time step"
  } else if (subtype == "peatland") {
    x <- magpie4::PeatlandArea(gdx, level = "cell", sum = FALSE)
    unit <- "Mha"
    description <- "peatland area by category"
  } else if (subtype == "landConservation") {
    # land conservation of the model, its policy countries and the fader of future conservation
    x <- list(landConservation = magpie4::landConservation(gdx, level = "cell"),
              policyCountries = gdx2::readGDX(gdx, "policy_countries22"),
              fader = gdx2::readGDX(gdx, "p22_conservation_fader", format = "first_found"))
    return(list(x = x, class = "list"))
  } else if (subtype == "snv") {
    # SNV share of each country (s29_snv_shr in the policy countries, s29_snv_shr_noselect in the
    # others), its fader, and the shares of the relocation targets; module 29 (cropland), or module 30 in
    # runs before the cropland module
    .first <- function(module29, module30) {
      return(gdx2::readGDX(gdx, module29, module30, format = "first_found"))
    }
    iso <- gdx2::readGDX(gdx, "iso")
    shareIso <- stats::setNames(rep(as.numeric(.first("s29_snv_shr_noselect", "s30_snv_shr_noselect")),
                                    length(iso)), iso)
    shareIso[.first("policy_countries29", "policy_countries30")] <- as.numeric(.first("s29_snv_shr", "s30_snv_shr"))
    x <- list(shareIso = shareIso,
              fader = .first("i29_snv_scenario_fader", "p30_snv_scenario_fader"),
              relocationX1 = as.numeric(.first("s29_snv_relocation_data_x1", "s30_snv_relocation_data_x1")),
              relocationX2 = as.numeric(.first("s29_snv_relocation_data_x2", "s30_snv_relocation_data_x2")))
    return(list(x = x, class = "list"))
  } else if (subtype == "bii") {
    x <- magpie4::BII(gdx, level = "cell", landClass = "all")
    unit <- "1"
    description <- "BII coefficients by land class"
  } else if (subtype == "climateFixYear") {
    return(list(x = as.numeric(gdx2::readGDX(gdx, "sm_fix_cc")), class = "numeric"))
  } else {
    stop("Unknown subtype '", subtype, "' in readMagpieFulldataGdx")
  }

  attr(x, "gdxMetadata") <- NULL

  return(list(x = x, min = 0, unit = unit, description = description))
}
