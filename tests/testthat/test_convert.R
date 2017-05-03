library('Meiosis')
context('convert')

data('exdat', package = 'Meiosis')

test_that("convert", {

ind <- structure(list(paternal = list(list(c(65L, 77L, 21L, 65L, 55L
), c(44.2292825428543, 79.8925440631938, 89.5366479233356, 120.316122102599,
157.515504024923)), list(c(77L, 55L, 21L, 65L, 21L, 55L, 77L,
21L, 65L), c(46.4656319277334, 92.0874250787024, 133.851070060911,
157.919032788446, 216.86555629969, 229.604085599934, 238.370587647335,
253.940346500535, 257.661027088761)), list(c(65L, 21L, 65L, 77L,
55L), c(1.77384141358994, 78.7934358375432, 99.8351353724615,
160.94833106676, 181.79538436234))), maternal = list(list(c(65L,
77L, 21L, 65L, 55L), c(23.8049165173269, 62.8291997013088, 89.5366479233356,
108.95375727535, 157.515504024923)), list(c(77L, 55L, 21L, 65L,
55L, 77L), c(46.4656319277334, 61.9024184249984, 133.851070060911,
142.240153051201, 229.604085599934, 257.661027088761)), list(
    c(65L, 21L, 77L, 21L, 65L, 77L, 55L), c(1.77384141358994,
    57.5905775935648, 69.9617374737411, 78.7934358375432, 155.335518385249,
    160.94833106676, 181.79538436234)))), .Names = c("paternal",
"maternal"))

ref <- list(list(c(1L, 1L, 0L, 1L, 0L, 0L, 1L, 0L, 1L), c(0L, 1L, 0L,
0L, 0L, 0L, 1L, 1L, 0L, 1L), c(1L, 0L, 0L, 0L, 1L)), list(c(1L,
1L, 1L, 1L, 0L, 0L, 1L, 0L, 1L), c(0L, 0L, 0L, 0L, 1L, 0L, 1L,
1L, 0L, 1L), c(1L, 0L, 0L, 0L, 1L)))

  co <- new(Converter, exdat$positions)
  alleles <- c(exdat$founder[[1L]][[1L]][[1L]], exdat$founder[[2L]][[1L]][[1L]])
  co$insert_founder(alleles, exdat$ind)
  alleles <- c(exdat$founder2[[1L]][[1L]][[1L]], exdat$founder2[[2L]][[1L]][[1L]])
  co$insert_founder(alleles, exdat$ind2)
  val <- co$convert(ind)
  testthat::expect_equal(val, expected = ref, tolerance = 1e-6)
})
