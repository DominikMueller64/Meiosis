#' Create founders with crossover-data.
#'
#' Create founders with crossover-data.
#' @param alleles Integer Vector of length 2, giving founder alleles.
#' @param L Numeric vector of chromosome lengths in cM.
#'
#' @return A nested list. The first level refers to gametes, the second
#' to chromosomes.
#'
#' @seealso \code{\link{cross_xo}}, \code{\link{meiosis_xo}}
#'
#' @author Dominik Mueller (\email{dominikmueller64@yahoo.de})
#'
#' @examples
#' create_xo_founder(alleles = c(-5L, 54L), L = c(32.2, 65.3, 88.2))
#'
#' @export
create_xo_founder <- function(alleles, L) {

  create_xo_gamete <- function(allele, L) {
    L <- as.numeric(L)
    allele <- as.integer(allele)
    if (!is.vector(L) || !is.numeric(L) || any(L <= 0.0))
      stop("'L' must be a numeric vector of (positive) lengths in cM.")

    if (!is.vector(allele) || !is.integer(allele) || length(allele) != 1L)
      stop("'allele' must be a single integer value.")

    lapply(X = L, FUN = function(l) list(allele, l))
  }

  if (length(alleles) != 2L)
    stop("'alleles' must have length 2.")

  lapply(X = alleles, FUN = function(allele) create_xo_gamete(allele, L))
}

