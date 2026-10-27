#' downscaleRun
#'
#' Downscale MAgPIE results from the given folder into that folder, as the tgz of the madrat collection
#' or, with archive = FALSE, as its files.
#'
#' @param outputfolder path to a folder containing fulldata.gdx and clustermap_*.rds,
#' the results are written here
#' @param model name of the full function, e.g. "DOWNSCALEDMAGPIE" or "ESM"
#' @param ... arguments passed on to retrieveData, e.g. rev, harmonizationPeriod, or with archive = FALSE to
#' the full function
#' @param cachefolder madrat cache folder for this call, e.g. in the output folder for cache files that only this
#' run reuses; NULL (default) keeps the madrat configuration
#' @param archive TRUE (default): retrieveData writes the collection as a tgz into the folder; FALSE: the files
#' of the collection go to the folder itself, without the tgz and the provenance files of retrieveData
#' (config.rds, renv.lock, diagnostics.log), and are computed again on each call
#' @param cores number of processes for the allocation kernel (option mrdownscale.cores, see
#' \code{\link{toolAllocateByCrossEntropy}}); above 1 not on Windows
#' @return Invisibly, the path to the newly created tgz archive, or the folder with archive = FALSE.
#'
#' @author Pascal Sauer, Patrick v. Jeetze
#' @export
downscaleRun <- function(outputfolder, model, ..., cachefolder = NULL, archive = TRUE, cores = 1) {
  outputfolder <- normalizePath(outputfolder)

  clustermap <- Sys.glob(file.path(outputfolder, "clustermap_*.rds"))
  gdx <- file.path(outputfolder, "fulldata.gdx")
  stopifnot(file.exists(gdx), length(clustermap) == 1)

  redirectSource("MagpieFulldataGdx", c(clustermap, gdx), linkOthers = FALSE)

  # the config and the 0.5 degree inputs of the run, without the downscaled outputs (cell.*), so that
  # the cache stays valid when outputs are added to the folder
  inputFiles <- Sys.glob(file.path(outputfolder, "*_0.5.mz"))
  runFiles <- c(file.path(outputfolder, "config.yml"), inputFiles[!startsWith(basename(inputFiles), "cell.")])
  redirectSource("MagpieRunFiles", runFiles[file.exists(runFiles)], linkOthers = FALSE)

  stopifnot(length(getConfig("redirections")) >= 1)
  localConfig(cachefolder = cachefolder)
  # an option, not an argument of the calc functions, so that the madrat cache does not depend on it
  withr::local_options(mrdownscale.cores = cores)
  if (archive) {
    return(retrieveData(model, ..., outputfolder = outputfolder, puc = FALSE))
  }
  # calcOutput(file = ) writes to the madrat output folder, the full function's other files to the working directory
  localConfig(outputfolder = outputfolder)
  fullFunction <- get(paste0("full", model))
  withr::with_dir(outputfolder, fullFunction(...))
  return(invisible(outputfolder))
}
