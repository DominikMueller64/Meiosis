.onAttach <- function(libname, pkgname) {
  packageStartupMessage(sprintf('Package Meiosis, version %s', packageVersion(pkgname)))
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

## Up until R 2.15.0, the require("methods") is needed but (now)
## triggers an warning from R CMD check
#.onLoad <- function(libname, pkgname){
#    #require("methods")  ## needed with R <= 2.15.0
#    loadRcppModules()
#}


## For R 2.15.1 and later this also works. Note that calling loadModule() triggers
## a load action, so this does not have to be placed in .onLoad() or evalqOnLoad().
loadModule("Module", TRUE)

