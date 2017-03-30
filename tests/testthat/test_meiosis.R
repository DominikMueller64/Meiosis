library('Meiosis')
context('meiosis')

data('exdat', package = 'Meiosis')
seed_rng(123L)

ref <-
structure(list(paternal = list(list(c(55L, 77L, 55L, 21L, 65L, 
21L, 77L, 55L), c(35.9549985909547, 46.2645149419177, 47.6562295838842, 
64.726958952326, 77.3588462323985, 89.5366479233356, 99.3885216073024, 
157.515504024923)), list(c(65L, 21L, 55L, 21L, 65L, 55L), c(113.002975706188, 
133.851070060911, 157.919032788446, 187.847696718557, 253.940346500535, 
257.661027088761)), list(c(55L, 21L, 65L, 21L, 77L), c(1.77384141358994, 
15.2210885765452, 77.7836075696138, 78.7934358375432, 181.79538436234
))), maternal = list(list(c(55L, 77L, 55L, 77L, 21L, 77L, 55L, 
21L, 65L), c(35.9549985909547, 46.2645149419177, 56.9873228330106, 
92.294090540975, 98.4319899227785, 99.3885216073024, 105.427297906101, 
122.866465852385, 157.515504024923)), list(c(77L, 55L, 65L, 21L, 
65L, 55L), c(44.8089322288614, 46.4656319277334, 113.002975706188, 
187.847696718557, 229.604085599934, 257.661027088761)), list(
    c(55L, 77L, 65L), c(77.5086992095846, 160.94833106676, 181.79538436234
    )))), .Names = c("paternal", "maternal"))


test_that("cross_xo", {
tmp <- cross_xo(exdat$founder, exdat$founder2, exdat$xoparam)
val <- cross_xo(tmp, tmp, exdat$xoparam)
testthat::expect_equal(val, expected = ref, tolerance = 1e-6)
})


seed_rng(123L)

ref <-
structure(list(paternal = list(c(0L, 0L, 0L, 1L, 0L, 1L, 0L, 
0L, 1L), c(0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L), c(0L, 0L, 
0L, 0L, 1L)), maternal = list(c(0L, 1L, 0L, 1L, 1L, 1L, 1L, 0L, 
0L), c(0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L), c(0L, 0L, 0L, 
0L, 1L))), .Names = c("paternal", "maternal"))

test_that("cross_geno", {
tmp <- cross_geno(exdat$ind, exdat$ind2, exdat$positions, exdat$xoparam)
val <- cross_geno(tmp, tmp, exdat$positions, exdat$xoparam)
testthat::expect_equal(val, expected = ref, tolerance = 1e-6)
})


