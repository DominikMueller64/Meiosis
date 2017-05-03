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

#' @name Converter
#' @docType class
#'
#' @title Converter reference class
#'
#' @description An object of type Converter is used to efficiently convert data from the
#' segmental to the genotypic representation. This class has two methods that should be
#' employed by the user: \code{insert_founder} and \code{convert}. Please look up their
#' documentation with \code{help(Converter$insert_founder)} and code{help(Converter$convert)}.
#' Before starting conversion, all founder alleles and founder genotypes have to be added via
#' the method \code{insert_founder}. Conversion is then done by the method \code{convert}.
#'
#' The constructor has two parameters, \code{positions} and \code{use_names}
#' (optional, defaults to \code{FALSE}).
#' Parameter \code{positions} must be a list of vectors with the genetic positions.
#' If these vectors are named, and if \code{use_names} is \code{TRUE},
#' the result of a conversion via \code{convert} will also be named.
#'
#' Please see the vignette (\code{vignette('Meiosis', package = 'Meiosis')}) for an example.
#'
NULL
## xyz <- setRefClass("xyz",
##             methods = list(
##               new = function(channel.labels=NA, green.channels=NA, avg.green.channel.label=NA) {
##                 "Create a new instance of the hvdvm class"
##                 NULL
##               },
##               digest = function(photo.conds.file.name, crop.files.path='./') {
##                 "Collect statistics from the raw images files sdfasdfsadf"
##                 NULL
##               }
##             ))

#' @name Converter$convert
#'
#' @title convert
#'
#' @description Convert an individual from its segmental to its genotypic representation.
#'
#' @param individual List. An individual to be converted.
NULL


#' @name Converter$insert_founder
#'
#' @title insert_founder
#'
#' @description Insert founder genotypes.
#'
#' @param keys Integer vector. The two (unique) founder alleles.
#' @param geno List. The genotypic data.
NULL


## ## destroy dummy class
## rm(xyz)
