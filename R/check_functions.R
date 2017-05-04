get_name <- function(x) {
  deparse(substitute(expr = x, env = parent.frame()))
}

is_sorted <- function(x) {
  !is.unsorted(x = x, strictly = TRUE)
}

check_ivec <- function(x, nm = NULL) {
  nm <- ifelse(test = is.null(nm), yes = get_name(x), no = nm)
  if (!is.integer(x) || !is.atomic(x))
    stop(sprintf('%s must be an integer vector.', nm))
}

check_locations <- function(x, nm = NULL) {
  nm <- ifelse(test = is.null(nm), yes = get_name(x), no = nm)
  if (!is.atomic(x) || !is.numeric(x) || !is_sorted(x) || x[1L] < 0.0) {
    msg <- sprintf(paste0('%s must be numeric, strictly increasingly sorted ',
                          'and contain only non-negative elements.'), nm)
    stop(msg)
  }
}

first_upper <- function(x) {
  paste0(toupper(substring(x, first = 1L, last = 1L)),
         substring(x, first = 2L, last = nchar(x)))
}

check_xo_gamete <- function(x, nm = NULL) {
  nm <- ifelse(test = is.null(nm), yes = get_name(x), no = nm)
  if (!is.list(x))
    stop(sprintf('Gamete %s must be itself a list.', nm))

  for (i in seq_along(x)) {
    core <- sprintf('chromosome %d on gamete %s', i, nm)
    a <- x[[i]]
    if (!is.list(a) || length(a) != 2L)
      stop(sprintf('%s must be itself a list of length 2.', first_upper(core)))
    check_ivec(x = a[[1L]], nm = sprintf('Founder alleles of %s', core))
    check_locations(x = a[[2L]], nm = sprintf('Locations of %s', core))
    if (length(a[[1L]]) != length(a[[2L]]))
      stop(sprintf('Founder alleles and locations of %s must have the same length.', core))
  }
}

check_geno_gamete <- function(x, nm = NULL) {
  nm <- ifelse(test = is.null(nm), yes = get_name(x), no = nm)
  if (!is.list(x))
    stop(sprintf('Gamete %s must be itself a list.', nm))
  for (i in seq_along(x)) {
    check_ivec(x = x[[i]], nm = sprintf('Alleles of chromosome %d on gamete %s', i, nm))
  }
}

#' @title Check genetic positions
#' @description Check if genetic positions are valid.
#'
#' @param x List. Genetic positions.
#'
#' @return \code{Null}
#' @export
check_positions <- function(x) {
  if (!is.list(x))
    stop(sprintf('%s must be a list.', get_name(x)))
  for (i in seq_along(x)) {
    check_locations(x[[i]], nm = sprintf('Positions of loci on chromosome %d', i))
  }
}

check_individual <- function(x, fun) {
  if (!is.list(x) || length(x) != 2L || length(unique(vapply(X = x, length, 1L))) != 1L)
    stop(sprintf('Individual %s must be a list of two gametes with the same length.',
                 get_name(x)))
  for (i in seq_along(x)) {
    g <- x[[i]]
    fun(x = g, nm = as.character(i))
  }
}

#' @title Check segmental representation of an individual.
#' @description Check if the segmental representation is valid.
#'
#' @param x List. Individual.
#'
#' @return \code{Null}
#'
#' @export
check_xo_individual <- function(x) {
  check_individual(x, fun = check_xo_gamete)
}

#' @title Check genotypic representation of an individual.
#' @description Check if the genotypic representation is valid.
#'
#' @param x List. Individual.
#'
#' @return \code{Null}
#'
#' @export
check_geno_individual <- function(x) {
  check_individual(x, fun = check_geno_gamete)
}

#' @title Check crossover parameters.
#' @description Check if crossover parameters are valid.
#'
#' @param x List. Crossover parameters.
#'
#' @return \code{Null}
#'
#' @export
check_xoparam <- function(x) {
  req <- c('L', 'm', 'p', 'obligate_chiasma', 'Lstar')
  if (!is.list(x) || !setequal(names(x), req)) {
    msg <- sprintf(paste0('Crossover parameters must be a list ',
                           'with elements %s.'), paste(req, collapse = ', '))
    stop(msg)
  }

  if (!is.atomic(x$L) || !is.double(x$L) || any(x$L <= 0.0))
    stop("'L' must be a numeric vector of (positive) chromosome lengths in CentiMorgan.")

  if (!is.integer(x$m) || length(x$m) != 1L || x$m < 0L)
    stop("'m' must be a non-negative integer.")

  if (!is.double(x$p) || length(x$p) != 1L || x$p < 0.0 || x$p > 1.0)
    stop("'p' must be a double between (including) 0 and 1")

  if (!is.logical(x$obligate_chiasma) || length(x$obligate_chiasma) != 1L)
    stop("'obligate_chiasma' must be a boolean value")

  if (!is.atomic(x$Lstar) || !is.double(x$Lstar) || any(x$Lstar <= 0.0 || any(x$Lstar > x$L)))
    stop("'Lstar' must be a numeric vector of (positive) chromosome lengths in CentiMorgan. ",
         "All its elements must be smaller than/equal to the corresponding elements in 'L'.")
}

