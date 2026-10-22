#' calcAgeClassesInput
#'
#' Age classes of an input dataset for the downscaling along its land transitions
#' (\code{\link{calcAgeClassesHighRes}}): unit stocks by land type and age class, and the age-specific flows
#' as operations in the order of the model.
#'
#' @details
#' Operations, each a list with \code{kind}:
#' \itemize{
#'   \item within: the unit flow \code{flowLr} of \code{pool} returns to the entry age classes, with the
#'   grid cell inflow of the land stages \code{inflowStages};
#'   \item transfer: the unit flow from \code{from} to \code{to} keeps its age; grid cell totals from the
#'   land stages \code{stages};
#'   \item shift: ageing by the length of the time step;
#'   \item exchange: the unit outflow \code{outflowLr} of \code{pool}, with grid cell totals from the land
#'   transitions outside \code{excludeStages}, the regrowth \code{regrowthLr} to the entry age classes, and
#'   the grid cell inflow to the entry age classes.
#' }
#' input = "magpie": the age classes of modules 32 and 35 (\code{magpie4::ageClasses},
#' \code{magpie4::ageClassTransitions}); per time step disturbance, recovery, ageing, maturation and the
#' outflow and regrowth of each land type, as in the presolve and the optimisation of the model.
#'
#' @param input name of an input dataset, currently only "magpie"
#' @return List with
#' \item{stocksLr}{magclass object of the unit stocks (Mha, names "land.type.ac")}
#' \item{pools}{land types with age classes, as labels "land.type" of \code{\link{calcLandTransitionsHighRes}}}
#' \item{groups}{named list of the output groups and their land types}
#' \item{classWidth}{years per age class}
#' \item{operations}{list of operations in model order}
#' \item{mapping}{data frame with columns cell and cluster}
#' @author Patrick v. Jeetze
#' @seealso \code{\link{calcAgeClassesHighRes}}
calcAgeClassesInput <- function(input) {
  if (input == "magpie") {
    stocksLr <- readSource("MagpieFulldataGdx", subtype = "ageClasses")
    flowsLr <- readSource("MagpieFulldataGdx", subtype = "ageClassTransitions")
    mapping <- readSource("MagpieFulldataGdx", subtype = "clustermap")

    # years per age class, from the names of the age classes (ac0, ac5, ...)
    ages <- as.numeric(sub("^ac", "", grep("^ac[0-9]+$", getItems(stocksLr, dim = "ac"), value = TRUE)))
    classWidth <- ages[2] - ages[1]

    # unit flow of a process and land type by age class
    .flow <- function(process, pool) {
      return(collapseDim(flowsLr[, , paste(process, pool, sep = ".")], dim = c(3.1, 3.2, 3.3)))
    }

    pools <- unique(paste(getItems(stocksLr, dim = 3.1, full = TRUE), getItems(stocksLr, dim = 3.2, full = TRUE),
                          sep = "."))
    forestry <- grep("^forestry[.]", pools, value = TRUE)
    # land types with disturbance (labels "land.type")
    disturbed <- grep("^disturbance[.]", getNames(flowsLr), value = TRUE)
    disturbed <- unique(sub("^disturbance[.](.*)[.][^.]+$", "\\1", disturbed))

    # presolve: disturbance, recovery, ageing, maturation; optimisation: outflow and regrowth
    operations <- list(list(kind = "within", pool = "secdforest.total",
                            flowLr = .flow("disturbance", "secdforest.total"), inflowStages = "disturbance"))
    for (pool in intersect(forestry, disturbed)) {
      operations[[length(operations) + 1]] <- list(kind = "within", pool = pool, flowLr = .flow("disturbance", pool))
    }
    operations[[length(operations) + 1]] <- list(kind = "transfer", from = "other.othernat", to = "other.youngsecdf",
                                                 flowLr = .flow("recovery", "other.othernat"),
                                                 stages = c("recoveryForestry", "recovery"))
    operations[[length(operations) + 1]] <- list(kind = "shift")
    operations[[length(operations) + 1]] <- list(kind = "transfer", from = "other.youngsecdf", to = "secdforest.total",
                                                 flowLr = .flow("maturation", "other.youngsecdf"),
                                                 stages = "maturation")
    for (pool in pools) {
      operations[[length(operations) + 1]] <- list(kind = "exchange", pool = pool,
                                                   outflowLr = .flow("outflow", pool),
                                                   regrowthLr = .flow("regrowth", pool),
                                                   excludeStages = c("disturbance", "recoveryForestry", "recovery",
                                                                     "maturation"))
    }

    out <- list(stocksLr = stocksLr, pools = pools,
                groups = list(secdforest = "secdforest.total", forestry = forestry,
                              other = c("other.othernat", "other.youngsecdf")),
                classWidth = classWidth, operations = operations, mapping = mapping)
  } else {
    stop("Unsupported input dataset \"", input, "\"")
  }

  return(list(x = out, class = "list", unit = "Mha",
              description = "Age classes of the input: unit stocks and age-specific flows as operations"))
}
