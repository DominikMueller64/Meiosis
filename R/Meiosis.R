#' Meiosis: Simulation of meiosis in plant breeding research
#'
#' This packages provides tools for the simulation of meiosis in plant breeding research.
#'
#' @docType package
#' @name Meiosis
NULL

#' @name seed_rng
#' @title Seed the random number generator.
#'
#' @description A Mersenne Twister pseudo-random generator of 32-bit numbers with a state size of
#' 10037 bits is used.
#'
#' @param seed Integer (default = \code{NULL}). If \code{NULL}, a random seed is used.
#'
#' @return Integer. The used seed.
#'
#' @examples
#' Meiosis::seed_rng(123L)
#' the_seed <- Meiosis::seed_rng()
#'
#' @export
NULL

#' @name calc_Lstar
#' @title Calculate adjusted chromosome length for obligate chiasma
#'
#' @description This function is an R-wrapper of an underlying C++ routine.
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
#' @examples
#' Meiosis::calc_Lstar(100, 0, 0)
#' Meiosis::calc_Lstar(60, 10, 0.1)
#'
#' @export
NULL

#' @name crossover
#' @title Simulate crossover locations.
#'
#' @description Simulate crossover locations on a single meiotic product using the
#' Stahl model.
#'
#' @param L Double. Length of the chromosome in cM. Must be > 50.
#' @param m Integer. Interference parameter for chi-square model. Must be non-negative.
#' @param p Double. Proportion of chiasmata coming from no-interference process
#' @param obligate_chiasma Logical. Should a least on the 4-strand bundle at meiosis.
#' @param Lstar Double. Reduced chromosome length as produced by \code{\link{calc_Lstar}}.
#'
#' @details Chiasma locations are a superposition of two
#' processes: a proportion p exhibiting no interference, and a
#' proportion (1-p) following the chi-square model with interference
#' parameter m.  Crossover locations are derived by thinning the
#' chiasma locations with probability 1/2.
#'
#' @details Simulations are under the Stahl model with the
#' interference parameter being an integer. This is an extension of
#' the chi-square model, but with chiasmata being the superposition of
#' two processes, one following the chi-square model and the other
#' exhibiting no interference.
#'
#' @details This function is an R-wrapper of an underlying C++ routine.
#' It is not intended for direct usage, but exposed for completeness.
#'
#' @author Dominik Mueller \email{dominikmueller64@@yahoo.de} and
#' Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#'
#' @return Numeric Vector. Crossover locations.
#'
#' @examples
#' Meiosis::crossover(300, 10, 0.5, FALSE, 300)
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
#' @export
NULL

#' @name meiosis_geno
#' @title Simulate meiosis (genotypic representation)
#'
#' @description Simulation of meiosis events.
#'
#' @param individual List. Individual.
#' @param positions List. Genetic positions.
#' @param xoparam List. Crossover parameters.
#' @param check Logical. Should checks be performed?
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
#' @description Simulation of a cross.
#'
#' @param father List. Father.
#' @param mother List. Mother.
#' @param positions List. Genetic positions.
#' @param xoparam List. Crossover parameters.
#' @param check Logical. Should checks be performed?
#'
#' @return List. A new individual.
#'
#' @examples
#' data(exdat, package = 'Meiosis')
#' Meiosis::cross_geno(exdat$ind, exdat$ind2, exdat$positions, exdat$xoparam)
#'
#' @export
NULL

#' @name self_geno
#' @title Produce selfing (genotypic representation)
#'
#' @description Simulation of a selfing.
#'
#' @param individual List. Individual.
#' @param positions List. Genetic positions.
#' @param xoparam List. Crossover parameters.
#' @param check Logical. Should checks be performed?
#'
#' @return List. A new individual.
#'
#' @examples
#' data(exdat, package = 'Meiosis')
#' Meiosis::self_geno(exdat$ind, exdat$positions, exdat$xoparam)
#'
#' @export
NULL


#' @name dh_geno
#' @title Produce doubled haploid (genotypic representation)
#'
#' @description Simulation of a doubled haploid.
#'
#' @param individual List. Individual.
#' @param positions List. Genetic positions.
#' @param xoparam List. Crossover parameters.
#' @param check Logical. Should checks be performed?
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
#' @description Simulation of meiosis events.
#'
#' @param individual List. Individual.
#' @param xoparam List. Crossover parameters.
#' @param check Logical. Should checks be performed?
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
#' @description Simulation of a cross.
#'
#' @param father List. Father.
#' @param mother List. Mother.
#' @param xoparam List. Crossover parameters.
#' @param check Logical. Should checks be performed?
#'
#' @return List. A new individual.
#'
#' @examples
#' data(exdat, package = 'Meiosis')
#' Meiosis::cross_xo(exdat$founder, exdat$founder2, exdat$xoparam)
#'
#' @export
NULL

#' @name self_xo
#' @title Produce selfing (segmental representation)
#'
#' @description Simulation of a selfing.
#'
#' @param individual List. Individual.
#' @param xoparam List. Crossover parameters.
#' @param check Logical. Should checks be performed?
#'
#' @return List. A new individual.
#'
#' @examples
#' data(exdat, package = 'Meiosis')
#' Meiosis::self_xo(exdat$founder, exdat$xoparam)
#'
#' @export
NULL


#' @name dh_xo
#' @title Produce doubled haploid (segmental representation)
#'
#' @description Simulation of a doubled haploid.
#'
#' @param individual List. Individual.
#' @param xoparam List. Crossover parameters.
#' @param check Logical. Should checks be performed?
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
#' @title Compute co-ancestry
#'
#' @description Compute the realized coefficient of co-ancestry between two individuals.
#'
#' @param individual_1 List. A first Individual.
#' @param individual_2 List. A second Individual. If \code{NULL}, self-relationship of
#' \code{individual_1} is computed.
#' @param check Logical. Should checks be performed?
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

#' @name realized_heter
#' @title Compute realized heterozygosity.
#'
#' @description Compute the realized heterozygosity of an individual, i.e., the proportion of
#' the genome where founder alleles match.
#'
#' @param individual List. Individual.
#' @param check Logical. Should checks be performed?
#'
#' @return double. Realized heterozygosity.
#'
#' @examples
#' data(exdat, package = 'Meiosis')
#' Meiosis::realized_heter(Meiosis::cross_xo(exdat$founder, exdat$founder, exdat$xoparam))
#'
#' @export
NULL

#' @name exdat
#' @title Example data
#'
#' @description Example data for illustrating the functions in Meiosis.
#'
#' @docType data
#' @format A list with 6 items
#' \describe{
#'   \item{xoparam}{List with crossover parameters}
#'   \item{positions}{Genetic positions of loci}
#'   \item{ind}{Genotypic data of an individual}
#'   \item{ind2}{Genotypic data of another individual}
#'   \item{founder}{Segmental data of an individual}
#'   \item{founder2}{Segmental data of another individual}
#' }
#' For details on the data structures, see \code{vignette('Meiosis', package = 'Meiosis')}.
"exdat"
