#' Meiosis: Simulation of meiosis in plant breeding research
#'
#' This packages provides tools for the simulation of meiosis in plant breeding research.
#'
#' @docType package
#' @name Meiosis
NULL

#' @name seed_rng
#' @title Seed the random number generator of the C++ routines.
#'
#' @description A Mersenne Twister pseudo-random generator of 32-bit numbers
#' with a state size of 10037 bits is used in the underlying C++ routines.
#'
#' @param seed Integer. If \code{NULL}, a random seed is used.
#'
#' @details The state of the random number generator used by the C++ routines is not
#' shared with the \code{R} session from which these routines are interfaced.
#'
#' @return Integer. The seed used.
#'
#' @examples
#' Meiosis::seed_rng(123L)
#' the_seed <- Meiosis::seed_rng()
#'
#' @export
NULL

#' @name calc_Lstar
#' @title Calculate adjusted chromosome length for an obligate chiasma.
#'
#' @description Calculate the reduced chromosome length that will give the target
#' expected number of chiasmata when conditioning on there being at
#' least one chiasma on the four-strand bundle.
#'
#' @inheritParams create_xoparam
#'
#' @param L Double. Length of the chromosome in cM. Must be > 50.
#' @param epsilon Double. The precision for finding the adjusted chromosome length.
#' Defaults to \code{NULL}, where a reasonably "high" precision is used.
#'
#' @details This function is an R-wrapper of an underlying C++ routine.
#' It is not intended for direct usage, but exposed for completeness.
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
#' @inheritParams create_xoparam
#'
#' @param L Double. Length of the chromosome in cM.
#' @param Lstar Double. Reduced chromosome length as produced by \code{\link{calc_Lstar}}.
#'
#' @inheritSection create_xoparam Model
#'
#' @details This function is an R-wrapper of an underlying C++ routine.
#' It is not intended for direct usage, but exposed for completeness.
#'
#' @return Double Vector. Crossover locations.
#'
#' @examples
#' Meiosis::crossover(300, 10, 0.5, FALSE, 300)
#'
#' @inherit create_xoparam references
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
#' @param use_names Logical. Should loci names be preserved if present?
#' Names in \code{individual} supersede the ones in \code{positions}.
#' @param check Logical. Should checks be performed?
#'
#' @return List. A new gamete.
#'
#' @examples
#' data(exdat, package = 'Meiosis')
#' Meiosis::meiosis_geno(exdat$ind, exdat$positions, exdat$xoparam, TRUE, TRUE)
#'
#' @export
meiosis_geno <- function(individual, positions, xoparam, use_names = FALSE, check = FALSE) {
  if (check) {
    check_geno_individual(individual)
    check_positions(positions)
    check_xoparam(xoparam)
  }
  .Call('Meiosis_meiosis_geno', PACKAGE = 'Meiosis', individual, positions, xoparam, use_names)
}


#' @name cross_geno
#' @title Cross individuals (genotypic representation)
#'
#' @description Simulation of a cross.
#'
#' @param father List. Father.
#' @param mother List. Mother.
#' @param positions List. Genetic positions.
#' @param xoparam List. Crossover parameters.
#' @param use_names Logical. Should loci names be preserved if present?
#' Names in \code{father} supersede the ones in \code{mother}, followed by \code{positions}.
#' @param check Logical. Should checks be performed?
#'
#' @return List. A new individual.
#'
#' @examples
#' data(exdat, package = 'Meiosis')
#' Meiosis::cross_geno(exdat$ind, exdat$ind2, exdat$positions, exdat$xoparam)
#'
#' @export
cross_geno <- function(father, mother, positions, xoparam, use_names = FALSE, check = FALSE) {
  if (check) {
    check_geno_individual(father)
    check_geno_individual(mother)
    check_positions(positions)
    check_xoparam(xoparam)
  }
    .Call('Meiosis_cross_geno', PACKAGE = 'Meiosis', father, mother, positions, xoparam, use_names)
}

#' @name self_geno
#' @title Produce selfing (genotypic representation)
#'
#' @description Simulation of a selfing.
#'
#' @param individual List. Individual.
#' @param positions List. Genetic positions.
#' @param xoparam List. Crossover parameters.
#' @param use_names Logical. Should loci names be preserved if present?
#' Names in \code{individual} supersede the ones in \code{positions}.
#' @param check Logical. Should checks be performed?
#'
#' @return List. A new individual.
#'
#' @examples
#' data(exdat, package = 'Meiosis')
#' Meiosis::self_geno(exdat$ind, exdat$positions, exdat$xoparam)
#'
#' @export
self_geno <- function(individual, positions, xoparam, use_names = FALSE, check = FALSE) {
  if (check) {
    check_geno_individual(individual)
    check_positions(positions)
    check_xoparam(xoparam)
  }
  .Call('Meiosis_self_geno', PACKAGE = 'Meiosis', individual, positions, xoparam, use_names)
}


#' @name dh_geno
#' @title Produce doubled haploid (genotypic representation)
#'
#' @description Simulation of a doubled haploid.
#'
#' @param individual List. Individual.
#' @param positions List. Genetic positions.
#' @param xoparam List. Crossover parameters.
#' @param use_names Logical. Should loci names be preserved if present?
#' Names in \code{individual} supersede the ones in \code{positions}.
#' @param check Logical. Should checks be performed?
#'
#' @return List. A new individual.
#'
#' @examples
#' data(exdat, package = 'Meiosis')
#' Meiosis::dh_geno(exdat$ind, exdat$positions, exdat$xoparam)
#'
#' @export
dh_geno <- function(individual, positions, xoparam, use_names = FALSE, check = FALSE) {
  if (check) {
    check_geno_individual(individual)
    check_positions(positions)
    check_xoparam(xoparam)
  }
    .Call('Meiosis_dh_geno', PACKAGE = 'Meiosis', individual, positions, xoparam, use_names)
}

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
meiosis_xo <- function(individual, xoparam, check = FALSE) {
  if (check) {
    check_xo_individual(individual)
    check_xoparam(xoparam)
  }
  .Call('Meiosis_meiosis_xo', PACKAGE = 'Meiosis', individual, xoparam)
}

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
cross_xo <- function(father, mother, xoparam, check = FALSE) {
  if (check) {
    check_xo_individual(father)
    check_xo_individual(mother)
    check_xoparam(xoparam)
  }
  .Call('Meiosis_cross_xo', PACKAGE = 'Meiosis', father, mother, xoparam)
}

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
self_xo <- function(individual, xoparam, check = FALSE) {
  if (check) {
    check_xo_individual(individual)
    check_xoparam(xoparam)
  }
  .Call('Meiosis_self_xo', PACKAGE = 'Meiosis', individual, xoparam)
}

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
dh_xo <- function(individual, xoparam, check = FALSE) {
  if (check) {
    check_xo_individual(individual)
    check_xoparam(xoparam)
  }
  .Call('Meiosis_dh_xo', PACKAGE = 'Meiosis', individual, xoparam)
}

#' @name realized_coancestry
#' @title Compute co-ancestry
#'
#' @description Compute the realized coefficient of co-ancestry between two individuals.
#'
#' @param individual_1 List. A first Individual in segmental representation.
#' @param individual_2 List. A second Individual in segmental representation.
#' If \code{NULL}, self-co-ancestry of \code{individual_1} is computed.
#' 
#' @param check Logical. Should checks be performed?
#'
#' @return Double. Realized coefficient of co-ancestry.
#'
#' @details The realized coefficient of co-ancestry is (herein) defined as the probability
#' that two alleles randomly drawn from one of the homologous chromosomes, each from one of the
#' individuals, are identical-by-descent (as specified by origins of founder alleles).
#'
#' @examples
#' \dontrun{
#' data(exdat, package = 'Meiosis')
#' Meiosis::realized_coancestry(Meiosis::dh_xo(exdat$founder, exdat$xoparam),
#'                              Meiosis::dh_xo(exdat$founder, exdat$xoparam))
#' }
#'
#' @export
realized_coancestry <- function(individual_1, individual_2 = NULL, check = FALSE) {
  if (check) {
    check_xo_individual(individual_1)
    check_xo_individual(individual_2)
  }
  .Call('Meiosis_realized_coancestry', PACKAGE = 'Meiosis', individual_1, individual_2)
}

#' @name realized_heter
#' @title Compute realized heterozygosity.
#'
#' @description Compute the realized heterozygosity of an individual, i.e., the proportion of
#' the genome where founder alleles match.
#'
#' @param individual List. Individual in segmental representation.
#' @param check Logical. Should checks be performed?
#'
#' @return double. Realized heterozygosity.
#'
#' @examples
#' \dontrun{
#' data(exdat, package = 'Meiosis')
#' Meiosis::realized_heter(Meiosis::cross_xo(exdat$founder, exdat$founder, exdat$xoparam))
#' }
#'
#' @export
realized_heter <- function(individual, check = FALSE) {
  if (check) {
    check_xo_individual(individual)
  }
  .Call('Meiosis_realized_heter', PACKAGE = 'Meiosis', individual)
}

#' @name exdat
#' @title Example data
#'
#' @description Example data for illustration.
#'
#' @docType data
#' @format A list with 6 items
#' \describe{
#'   \item{xoparam}{List with crossover parameters (see \code{\link{create_xoparam}})}
#'   \item{positions}{Genetic positions of loci}
#'   \item{ind}{Genotypic data of an individual}
#'   \item{ind2}{Genotypic data of another individual}
#'   \item{founder}{Segmental data of an individual}
#'   \item{founder2}{Segmental data of another individual}
#' }
#' For details on the data structures, see \code{vignette('Introduction', package = 'Meiosis')}.
"exdat"
