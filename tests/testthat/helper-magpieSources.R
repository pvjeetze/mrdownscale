# Replaces readSource and calcOutput in the calc functions of the MAgPIE downscaling: readSource
# returns sources[["<type>.<subtype>"]] or else sources[[subtype]]; calcOutput returns outputs[[key]],
# with key the calc name and its arguments other than NULL and the defaults of the calc function (for
# example "LandTransitionsHighRes.input=magpie"), or else runs the calc function.
localMagpieSources <- function(sources, outputs = list(), env = parent.frame()) {
  local_mocked_bindings(
    readSource = function(type, subtype = NULL, ...) {
      for (key in c(paste0(type, ".", subtype), subtype)) {
        if (key %in% names(sources)) {
          return(sources[[key]])
        }
      }
      stop("no synthetic source for ", type, " subtype ", subtype)
    },
    calcOutput = function(type, ..., aggregate = TRUE) {
      args <- list(...)
      calc <- get(paste0("calc", type))
      defaults <- formals(calc)
      keyArgs <- character(0)
      for (name in names(args)) {
        if (!is.null(args[[name]]) && !identical(args[[name]], defaults[[name]])) {
          keyArgs <- c(keyArgs, paste0(name, "=", paste(args[[name]], collapse = ",")))
        }
      }
      key <- paste(c(type, keyArgs), collapse = ".")
      if (key %in% names(outputs)) {
        return(outputs[[key]])
      }
      return(do.call(calc, args)[["x"]])
    },
    .env = env
  )
}
