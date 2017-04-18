#' @rdname seed_rng
seed_rng <- function(seed = NULL) {
    .Call('Meiosis_seed_rng', PACKAGE = 'Meiosis', seed)
}

#' @rdname calc_Lstar
calc_Lstar <- function(L, m, p, epsilon = NULL) {
    .Call('Meiosis_calc_Lstar', PACKAGE = 'Meiosis', L, m, p, epsilon)
}

#' @rdname crossover
crossover <- function(L, m, p, obligate_chiasma, Lstar) {
    .Call('Meiosis_crossover', PACKAGE = 'Meiosis', L, m, p, obligate_chiasma, Lstar)
}

#' @rdname meiosis_geno
meiosis_geno <- function(individual, positions, xoparam, check = FALSE) {
  if (check) {
    check_geno_individual(individual)
    check_positions(positions)
    check_xoparam(xoparam)
  }
  .Call('Meiosis_meiosis_geno', PACKAGE = 'Meiosis', individual, positions, xoparam)
}

#' @rdname meiosis_xo
meiosis_xo <- function(individual, xoparam, check = FALSE) {
  if (check) {
    check_xo_individual(individual)
    check_xoparam(xoparam)
  }
  .Call('Meiosis_meiosis_xo', PACKAGE = 'Meiosis', individual, xoparam)
}

#' @rdname cross_geno
cross_geno <- function(father, mother, positions, xoparam, check = FALSE) {
  if (check) {
    check_geno_individual(father)
    check_geno_individual(mother)
    check_positions(positions)
    check_xoparam(xoparam)
  }
  .Call('Meiosis_cross_geno', PACKAGE = 'Meiosis', father, mother, positions, xoparam)
}

#' @rdname cross_xo
cross_xo <- function(father, mother, xoparam, check = FALSE) {
  if (check) {
    check_xo_individual(father)
    check_xo_individual(mother)
    check_xoparam(xoparam)
  }
  .Call('Meiosis_cross_xo', PACKAGE = 'Meiosis', father, mother, xoparam)
}

#' @rdname self_geno
self_geno <- function(individual, positions, xoparam, check = FALSE) {
  if (check) {
    check_geno_individual(individual)
    check_positions(positions)
    check_xoparam(xoparam)
  }
  .Call('Meiosis_self_geno', PACKAGE = 'Meiosis', individual, positions, xoparam)
}

#' @rdname self_xo
self_xo <- function(individual, xoparam, check = FALSE) {
  if (check) {
    check_xo_individual(individual)
    check_xoparam(xoparam)
  }
  .Call('Meiosis_self_xo', PACKAGE = 'Meiosis', individual, xoparam)
}

#' @rdname dh_geno
dh_geno <- function(individual, positions, xoparam, check = FALSE) {
  if (check) {
    check_geno_individual(individual)
    check_positions(positions)
    check_xoparam(xoparam)
  }
  .Call('Meiosis_dh_geno', PACKAGE = 'Meiosis', individual, positions, xoparam)
}

#' @rdname dh_xo
dh_xo <- function(individual, xoparam, check = FALSE) {
  if (check) {
    check_xo_individual(individual)
    check_xoparam(xoparam)
  }
  .Call('Meiosis_dh_xo', PACKAGE = 'Meiosis', individual, xoparam)
}

#' @rdname realized_heter
realized_heter <- function(individual, check = FALSE) {
  if (check) {
    check_xo_individual(individual)
  }
  .Call('Meiosis_realized_heter', PACKAGE = 'Meiosis', individual)
}

#' @rdname realized_coancestry
realized_coancestry <- function(individual_1, individual_2 = NULL, check = FALSE) {
  if (check) {
    check_xo_individual(individual_1)
    check_xo_individual(individual_2)
  }
  .Call('Meiosis_realized_coancestry', PACKAGE = 'Meiosis', individual_1, individual_2)
}
