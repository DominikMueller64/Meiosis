#' Create a parameter list for simulating crossovers.
#'
#'
#' Create a parameter list for simulating crossover locations on a single meiotic product
#' using the Stahl model.
#'
#' @param L Double vector. Chromosome lengths in CentiMorgan (cM).
#' @param m Integer. The interference parameter (\code{m = 0} is no interference).
#' @param p Double. Proportion of chiasmata from no-interference mechanism.
#' (\code{p = 0} gives pure chi-square model)
#' @param obligate_chiasma Logical. If TRUE, require an obligate chiasma on the
#' 4-strand bundle at meiosis. Only possible if all chromosomes are longer than 50 cM.
#'
#' @return A nested list as long as the number of chromosomes with parameters used for
#' simulating crossover events.
#'
#' @section Model:
#' Chiasma locations are a superposition of two
#' processes: a proportion p exhibiting no interference, and a
#' proportion \code{(1 - p)} following the chi-square model with interference
#' parameter m.  Crossover locations are derived by thinning the
#' chiasma locations with probability \code{1/2}.
#'
#' Simulations are under the Stahl model with the
#' interference parameter being an integer. This is an extension of
#' the chi-square model, but with chiasmata being the superposition of
#' two processes, one following the chi-square model and the other
#' exhibiting no interference.
#'
#' @details
#' If \code{obligate_chiasma = TRUE} and all chromsome are longer than 50 cM, reduced
#' chromosome length with be internally calculated (\code{calc_Lstar}) that will give the
#' target expected number of chiasmata when conditioning on there being at least one
#' chiasma on the four-strand bundle.
#'
#' This description is adopted from the \href{https://github.com/kbroman/simcross}{simcross}
#' package.
#'
#' @references
#' Copenhaver, G. P., Housworth, E. A. and Stahl, F. W. (2002) Crossover
#' interference in arabidopsis.  \emph{Genetics} \bold{160}, 1631--1639.
#'
#' Foss, E., Lande, R., Stahl, F. W. and Steinberg, C. M. (1993) Chiasma
#' interference as a function of genetic distance. \emph{Genetics}
#' \bold{133}, 681--691.
#'
#' Zhao, H., Speed, T. P. and McPeek, M. S. (1995) Statistical analysis
#' of crossover interference using the chi-square model.  \emph{Genetics}
#' \bold{139}, 1045--1056.
#'
#' @examples
#' create_xoparam(L = c(50.5, 100.1, 133.5))
#'
#' @export
create_xoparam <- function(L, m = 0L, p = 1.0,
                            obligate_chiasma = FALSE) {

  if (!is.numeric(p) || length(p) != 1L || p < 0.0 || p > 1.0)
    stop("'p' must be a double between (including) 0 and 1")

  if (!is.numeric(m) || length(p) != 1L || m < 0.0 || !isTRUE(all.equal(m, round(m))))
    stop("'m' must be an integer greater or equal to 0.")

  if (!is.atomic(L) || !is.double(L) || any(L <= 0.0))
    stop("'L' must be a numeric vector of (positive) chromosome lengths in CentiMorgan.")

  p <- as.double(p)
  m <- as.integer(round(m))
  L <- as.double(L)
  obligate_chiasma <- as.logical(obligate_chiasma)

  if (obligate_chiasma && any(L <= 50.0))
    stop("'obligate_chiasma = TRUE' is only possible if all chromosomes are longer than 50 cM.")

  ret <- list()
  ret$L <- L
  ret$m <- m
  ret$p <- p
  ret$obligate_chiasma <- obligate_chiasma
  if (obligate_chiasma)
    Lstar <- vapply(X = L, FUN = Meiosis::calc_Lstar, FUN.VALUE = numeric(1L),
                    m = m, p = p, epsilon = sqrt(.Machine$double.eps))
  else {
    Lstar <- L
  }
  ret$Lstar <- Lstar
  ret
}
