library('Meiosis')
context('convert')

data('exdat', package = 'Meiosis')
seed_rng(123L)

ref <-
list(list(c(0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 1L), c(0L, 1L, 0L, 
0L, 0L, 0L, 0L, 1L, 1L, 1L), c(0L, 0L, 0L, 0L, 1L)), list(c(0L, 
1L, 0L, 1L, 1L, 1L, 1L, 0L, 0L), c(0L, 1L, 0L, 0L, 0L, 0L, 0L, 
0L, 0L, 1L), c(0L, 0L, 0L, 0L, 1L)))

test_that("convert", {
co <- new(Converter, exdat$positions)
alleles <- c(exdat$founder[[1L]][[1L]][[1L]], exdat$founder[[2L]][[1L]][[1L]])
co$insert_founder(alleles, exdat$ind)
alleles <- c(exdat$founder2[[1L]][[1L]][[1L]], exdat$founder2[[2L]][[1L]][[1L]])
co$insert_founder(alleles, exdat$ind2)
tmp <- cross_xo(exdat$founder, exdat$founder2, exdat$xoparam)
val <- co$convert(cross_xo(tmp, tmp, exdat$xoparam))
testthat::expect_equal(val, expected = ref, tolerance = 1e-6)
})
