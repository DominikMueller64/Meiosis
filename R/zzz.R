.onAttach <- function(libname, pkgname) {
  packageStartupMessage(sprintf('Package Meiosis, version %s', utils::packageVersion(pkgname)))
  packageStartupMessage('Dominik Mueller')
  packageStartupMessage("Use Meiosis::seed_rng to seed the random number generator (Mersenne Twister, std::mt19937)")
  packageStartupMessage("Otherwise a random seed is used.")
  invisible()
}

.onLoad <- function(libname, pkgname) {
  seed_rng()
  invisible()
}

.onUnload <- function (libpath) {
  library.dynam.unload("Meiosis", libpath)
}

Rcpp::loadModule("Module", TRUE)

## Rcpp::setRcppClass("Converter",
##                    module = "Module",
##                    fields = list(more = "character"),
##                    methods = list(
##                      test = function(what) message("Testing: ", what, "; ", more)),
##                    saveAs = "genConverter"
##                    )
