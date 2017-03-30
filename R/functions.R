check_sorted <- function(x) {
  !is.unsorted(x = x, strictly = TRUE)
}

check_alleles <- function(x) {
  if (!is.integer(x) || !is.atomic(x)) stop("Alleles must be an integer vector.")
  return(TRUE)
}

check_locations <- function(x) {
  if (!check_sorted(x) || x[1L] < 0.0 || !is.numeric(x) || !is.atomic(x))
    stop(paste0("Locations must be numeric, increasingly sorted (strictly) ",
                "and contain only non-negative elements."))
  return(TRUE)
}

check_xo_gamete <- function(x) {
  if (!is.list(x)) stop("'gamete' must be a list.")
  for (a in x) {
    if (!is.list(a)) stop("Each element of 'gamete' must be itself a list.")
    check_alleles(a[[1L]])
    check_locations(a[[2L]])
  }
  return(TRUE)
}

check_geno_gamete <- function(x) {
  if (!is.list(x)) stop("'gamete' must be a list.")
  for (a in x) check_alleles(a)
  return(TRUE)
}

#' @title Check genetic positions
#' @description Check if genetic positions are valid.
#'
#' @param x List. Genetic positions.
#' @export
check_positions <- function(x) {
  if (!is.list(x)) stop("'x' must be a list.")
  for (p in x) check_locations(p)
  return(TRUE)
}

#' @title Check segmental representation.
#' @description Check if the segmental representation is valid.
#'
#' @param x List. Individual
#' @export
check_xo_individual <- function(x) {
  if (!is.list(x) || length(x) != 2L)
    stop("An individuals must be a 'list' of length two (two gametes).")
  for (g in x) check_xo_gamete(g)
  return(TRUE)
}

#' @title Check genotypic representation.
#' @description Check if the genotypic representation is valid.
#'
#' @param x List. Individual
#' @export
check_geno_individual <- function(x) {
  if (!is.list(x) || length(x) != 2L)
    stop("An individuals must be a 'list' of length two (two gametes).")
  for (g in x) check_geno_gamete(g)
  return(TRUE)
}


