#' Meiosis: A package for the simulation of meiosis events in diploids.
#'
#' This packages provides tools for the simulation of meiosis events in diploid species.
#'
#' @docType package
#' @name Meiosis
NULL

#' @name seed_rng
#' @title Seed the random number generator.
#'
#' A Mersenne Twister pseudo-random generator of 32-bit numbers with a state size of
#' 10037 bits is used.
#'
#' @param seed Integer (default = \code{NULL}). If \code{NULL}, a random seed is used.
#'
#'
#' @return Integer. The used seed.
#'
#' @seealso \url{http://www.cplusplus.com/reference/random/mt19937/}
#'
#' @examples
#' seed_rng(123L)
#' the_seed <- seed_rng()
#'
#' @export
NULL

#' @name calc_Lstar
#' @title Calculate adjusted chromosome length for obligate chiasma
#'
#' @details This function is an R-wrapper of an underlying C++ routine.
#' It is not intended for direct usage, but exposed for completeness.
#'
#' Calculate the reduced chromosome length that will give the target
#' expected number of chiasmata when conditioning on there being at
#' least one chiasma on the four-strand bundle.
#'
#' @param L Double. Length of the chromosome in cM. Must be > 50.
#' @param m Integer. Interference parameter for chi-square model. Must be non-negative.
#' @param p Double. Proportion of chiasmata coming from no-interference process
#' @param epsilon Double. The precision for finding the adjusted chromosome length.
#' Defaults to \code{NULL}, where a "high" precision is used.
#'
#' @return Double. Adjusted chromosome length.
#'
#'
#' @author Orginal R function: Karl Broman, C++ routine: Dominik Mueller.
#'
#' @seealso \code{\link[simcross]{calc_Lstar}}
#'
#' @examples
#' calc_Lstar(100, 0, 0)
#' calc_Lstar(60, 10, 0.1)
#'
#' @export
NULL

#' @name crossover
#' @title Simulate crossover locations.
#'
#' @details This function is an R-wrapper of an underlying C++ routine.
#' It is not intended for direct usage, but exposed for completeness.
#'
#' @param L Double. Length of the chromosome in cM. Must be > 50.
#' @param m Integer. Interference parameter for chi-square model. Must be non-negative.
#' @param p Double. Proportion of chiasmata coming from no-interference process
#' @param obligate_chiasma Logical. Should a least on the 4-strand bundle at meiosis.
#' @param Lstar Double. Reduced chromosome length as produced by \code{\link{calc_Lstar}}.
#'
#' @return Numeric Vector. Crossover locations.
#'
#' @seealso \code{\link[simcross]{sim_crossover}}
#'
#' @examples
#' crossover(300, 10, 0.5, FALSE, 300)
#'
#' @export
NULL

#' @name meiosis_geno
#' @title Simulate meiosis (genotypic representation)
#'
#' @param individual List. Individual.
#' @param positions List. Genetic positions.
#' @param xoparam List. Crossover parameters.
#'
#' @return List. A new gamete.
#'
#' @examples
#' data(exdat, package = 'Meiosis')
#' Meiosis::meiosis_geno(exdat$ind, exdat$positions, exdat$xoparam)
#'
#' @export
NULL


#' @name cross_geno
#' @title Cross individuals (genotypic representation)
#'
#' @param father List. Father.
#' @param mother List. Father.
#' @param positions List. Genetic positions.
#' @param xoparam List. Crossover parameters.
#'
#' @return List. A new individual.
#'
#' @examples
#' data(exdat, package = 'Meiosis')
#' Meiosis::cross_geno(exdat$ind, exdat$ind2, exdat$positions, exdat$xoparam)
#'
#' @export
NULL


#' @name dh_geno
#' @title Produce doubled haploid (genotypic representation)
#'
#' @param individual List. Individual.
#' @param positions List. Genetic positions.
#' @param xoparam List. Crossover parameters.
#'
#' @return List. A new individual.
#'
#' @examples
#' data(exdat, package = 'Meiosis')
#' Meiosis::dh_geno(exdat$ind, exdat$positions, exdat$xoparam)
#'
#' @export
NULL

#' @name meiosis_xo
#' @title Simulate meiosis (segmental representation)
#'
#' @param individual List. Individual.
#' @param xoparam List. Crossover parameters.
#'
#' @return List. A new gamete.
#'
#' @examples
#' data(exdat, package = 'Meiosis')
#' Meiosis::meiosis_xo(exdat$founder, exdat$xoparam)
#'
#' @export
NULL


#' @name cross_xo
#' @title Cross individuals (segmental representation)
#'
#' @param father List. Father.
#' @param mother List. Father.
#' @param xoparam List. Crossover parameters.
#'
#' @return List. A new individual.
#'
#' @examples
#' data(exdat, package = 'Meiosis')
#' Meiosis::cross_xo(exdat$founder, exdat$founder2, exdat$xoparam)
#'
#' @export
NULL


#' @name dh_xo
#' @title Produce doubled haploid (segmental representation)
#'
#' @param individual List. Individual.
#' @param positions List. Genetic positions.
#' @param xoparam List. Crossover parameters.
#'
#' @return List. A new individual.
#'
#' @examples
#' data(exdat, package = 'Meiosis')
#' Meiosis::dh_xo(exdat$founder, exdat$xoparam)
#'
#' @export
NULL

#' @name realized_coancestry
#' @title Compute the realized coefficient of co-ancestry between two individuals.
#'
#' @param individual_1 List. A first Individual.
#' @param individual_2 List. A second Individual. If \code{NULL}, self-relationship of
#' \code{individual_1} is computed.
#'
#' @return double. Realized coefficient of co-ancestry.
#'
#' @details The realized coefficient of co-ancestry is (herein) defined as the probability
#' that two alleles randomly drawn from one of the homologous chromosomes, each from one of the
#' individuals, are identical-by-descent (as specified by origins of founder alleles).
#'
#' @examples
#' data(exdat, package = 'Meiosis')
#' Meiosis::realized_coancestry(Meiosis::dh_xo(exdat$founder, exdat$xoparam),
#'                              Meiosis::dh_xo(exdat$founder, exdat$xoparam))
#'
#' @export
NULL

#' @name exdat
#' @docType data
#' Example data
#'
#' Example data for illustrating the functions in Meiosis.
#'
#' @format A list with 6 items
#' \describe{
#'   \item{xoparam}{List with crossover parameters}
#'   \item{positions}{Genetic positions of loci}
#'   \item{ind}{Genotypic data of an individual}
#'   \item{ind2}{Genotypic data of another individual}
#'   \item{founder}{Segmental data of an individual}
#'   \item{founder2}{Segmental data of another individual}
#' }
#'
#' For details on the data structures, see \code{vignette('Meiosis', package = 'Meiosis')}.
"exdat"
