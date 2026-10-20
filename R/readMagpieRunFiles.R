#' readMagpieRunFiles
#'
#' Read function for the config and the 0.5 degree input files that MAgPIE exports to the output
#' folder of a run (see \code{\link{downscaleRun}}).
#'
#' @param subtype type of data to be read in: config, landuseInit, availableCropland, potentialForest,
#' potentialForestUncorrected, urban, wdpa, conservationPriority, snvTarget, peatland, sideLayers, or
#' available (the subtypes whose file exists)
#' @author Patrick v. Jeetze
readMagpieRunFiles <- function(subtype) {
  file <- c(config = "config.yml",
            landuseInit = "avl_land_full_t_0.5.mz",
            availableCropland = "avl_cropland_0.5.mz",
            potentialForest = "pot_forest_area_0.5.mz",
            potentialForestUncorrected = "pot_forest_area_uncorrected_0.5.mz",
            urban = "f34_urbanland_0.5.mz",
            wdpa = "wdpa_baseline_0.5.mz",
            conservationPriority = "consv_prio_areas_0.5.mz",
            snvTarget = "SNVTargetCropland_0.5.mz",
            peatland = "f58_peatland_area_0.5.mz",
            sideLayers = "luh2_side_layers_0.5.mz")

  if (subtype == "available") {
    return(list(x = names(file)[file.exists(file)], class = "character"))
  }
  if (!subtype %in% names(file)) {
    stop("Unknown subtype '", subtype, "' in readMagpieRunFiles")
  }
  if (!file.exists(file[[subtype]])) {
    stop(file[[subtype]], " not found in the output folder of the run")
  }
  if (subtype == "config") {
    return(list(x = gms::loadConfig(file[["config"]]), class = "list"))
  }

  unit <- if (subtype == "sideLayers") "1" else "Mha"
  return(list(x = read.magpie(file[[subtype]]), unit = unit, description = paste("MAgPIE run input", file[[subtype]])))
}
