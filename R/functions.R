check_sorted <- function(x) {
  !is.unsorted(x = x, strictly = TRUE)
}

check_alleles <- function(x) {
  if (!is.integer(x) || !is.atomic(x)) stop("Alleles must be an integer vector.")
}

check_locations <- function(x) {
  if (!check_sorted(x) || x[1L] < 0.0 || !is.numeric(x) || !is.atomic(x))
    stop(paste0("Locations must be numeric, increasingly sorted (strictly) ",
                "and contain only non-negative elements."))
}

#' @export
check_pos <- function(x) {
  if (!is.list(x)) stop("'pos' must be a list.")
  for (p in x) check_locations(p)
}

#' @export
check_xodat_gamete <- function(x) {
  if (!is.list(x)) stop("'gamete' must be a list.")
  for (a in x) {
    if (!is.list(a)) stop("Each element of 'gamete' must be itself a list.")
    check_alleles(a[[1L]])
    check_locations(a[[2L]])
  }
}

#' @export
check_xodat_individual <- function(x) {
  if (!is.list(x) || length(x) != 2L)
    stop("An individuals must be a 'list' of length two (two gametes).")
  for (g in x) check_xodat_gamete(g)
}
