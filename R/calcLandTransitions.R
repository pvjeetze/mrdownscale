#' calcLandTransitions
#'
#' Ex-Post estimation of land use transitions based on land use state
#' information
#'
#' @param outputFormat format in which the outputs should be prepared. Currently,
#' only "ESM" for earth system model compatible input data is available.
#' @inheritParams calcNonlandInput
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @param yearsSubset vector of years to keep in the output dataset
#' @param harmonization name of harmonization method, see \code{\link{toolGetHarmonizer}}
#' @param downscaling name of downscaling method: "proportional" (\code{luscale::interpolate2}) or "transitions"
#' (\code{\link{toolDownscaleTransitions}})
#' @param gross either boolean or a magpie object containing bidirectional
#' transition shares relative to the area of the involved land pools (transition
#' divided by the area of the land pool in the "from" sub dimension). If set to
#' FALSE only net transitions will be returned. If set to TRUE an internal
#' gross transition estimate based on average gross transitions in LUH2 in the
#' period from 1995 to 2015 will be used.
#' @return land use transition data
#' @author Jan Philipp Dietrich, Pascal Sauer
calcLandTransitions <- function(outputFormat, input, harmonizationPeriod, yearsSubset,
                                harmonization, downscaling, gross) {
  if (outputFormat != "ESM") {
    stop("Can only report for outputFormat = 'ESM'")
  }

  land <- calcOutput("LandReport", outputFormat = "ESM", input = input,
                     harmonizationPeriod = harmonizationPeriod, yearsSubset = yearsSubset,
                     harmonization = harmonization, downscaling = downscaling, aggregate = FALSE)
  land <- land[, , grep("(_|manaf)", getItems(land, dim = 3), invert = TRUE, value = TRUE)]

  # add extra year as copy of last year to get gross transitions (net zero) for
  # the time after the covered period
  lastYear <- max(getYears(land, as.integer = TRUE))
  land <- mbind(land, setYears(land[, lastYear, ], lastYear + 1))

  timestepLengths <- new.magpie(years = getYears(land)[-nyears(land)])
  timestepLengths[, , ] <- diff(getYears(land, as.integer = TRUE))

  l <- 2
  sequence <- seq(1, nyears(land) - 1, l)
  tempfolder <- withr::local_tempdir()
  for (i in sequence) {
    i2 <- min(nyears(land), i + l)
    transition <- toolMinimalTransitions(land[, i:i2, ], gross = gross)

    # check if the calculated transition is consistent with the states (land)
    netChangeLand <- setYears(land[, (i + 1):i2, ], getYears(land)[i:(i2 - 1)]) - land[, i:(i2 - 1), ]

    netChangeTrans <- -dimSums(transition, dim = "to")
    to <- getItems(transition, "to")
    netChangeTrans[, , to] <- netChangeTrans[, , to] + dimSums(transition, "from")
    getYears(netChangeTrans) <- getYears(netChangeLand)
    netChangeTrans <- netChangeTrans * timestepLengths[, getYears(netChangeTrans), ]

    if (max(transition) > 1.0001) {
      toolStatusMessage("warn", paste0("Some shares are > 1 (max: ", max(transition), "), setting those to 1"))
    }
    transition[transition > 1] <- 1

    toolExpectLessDiff(netChangeTrans, netChangeLand[, , getItems(netChangeTrans, 3)],
                       10^-6, paste0("Transitions are consistent to states from ",
                                     getYears(land)[i], " to ", getYears(land)[i2]))

    write.magpie(transition, paste0(tempfolder, "/", i, ".mz"))
  }

  # allocate enough memory for everything instead of mbinding to reduce memory requirement
  message("allocating memory...")
  coords <- getItems(land, dim = 1)
  years <- utils::head(getYears(land, as.integer = TRUE), -1) + 1
  rm(land)
  variables <- getItems(read.magpie(paste0(tempfolder, "/", sequence[1], ".mz")), dim = 3)
  out <- list(x = new.magpie(coords, years, variables, fill = NA_real_,
                             sets = c("x", "y", "year", "from", "to")),
              isocountries = FALSE,
              unit = "1",
              min = 0,
              max = 1,
              description = paste("MAgPIE land use transition data estimated from downscaled land use ",
                                  "state information. Unit is share of cell area."))

  for (j in seq_along(sequence)) {
    i <- sequence[j]
    message(j, "/", length(sequence), " combining...")
    y <- years[i:min(length(years), i + l - 1)]
    out$x[, y, ] <- read.magpie(paste0(tempfolder, "/", i, ".mz"))
  }
  message("combining done, returning")

  return(out)
}
