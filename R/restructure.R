#' @title Restructure genotypic data
#'
#' @description Restructure genotypic data into several forms.
#'
#' This function is not yet official part of the package:
#' @keywords internal
#'
#' @param individual An individual (genotypic representation).
#' @param struc Character. Either of \code{c('_', 'g', 'c', 'cg')}
#' @param check Logical. Should checks be performed?
#'
#' @details If \code{struc == '_'} (default), the result is a matrix (both chromosomes
#' and gametes are concatenated). If \code{struc == 'g'}, only gametes are concatenated
#' (chromosomes remain separated). If \code{struc == 'c'}, only chromosomes are concatenated
#' (gametes remain separated), and if \code{struc == 'cg'}, gametes get nested
#' within chromosomes, everything being separated.
#'
#' @return A matrix or a list, depending on the parameter \code{struc}.
#'
#' @examples
#' \dontrun{
#' data('exdat', package = 'Meiosis')
#' str(restructure(exdat$ind))
#' str(restructure(exdat$ind, struc = 'c'))
#' str(restructure(exdat$ind, struc = 'g'))
#' str(restructure(exdat$ind, struc = 'cg'))
#' }
#'
restructure <- function(individual, struc = c('_', 'g', 'c', 'cg'),
                        check = FALSE) {

  struc <- match.arg(arg = struc)

  if (check) {
    check_geno_individual(x = individual)
  }

  if (struc %in% c('_', 'g')) {
    ## Faster than use.names = TRUE, no checking for duplicate names.
    ind <- lapply(X = individual, FUN = function(x) {
      nm <- unlist(x = lapply(X = x, FUN = names), recursive = FALSE,
                   use.names = FALSE)
      object <- unlist(x, recursive = FALSE, use.names = FALSE)
      names(object) <- nm
      object
    })

    if (struc == '_') ind <- do.call(what = rbind, args = ind)
  }
  else if (struc %in% c('c', 'cg')) {
    ind <- lapply(X = seq_along(individual[[1L]]),
                  FUN = function(i) lapply(X = individual, FUN = `[[`, i = i))

    if (struc == 'c') ind <- lapply(X = ind, FUN = do.call, what = rbind)
  }
  return(ind)
}
