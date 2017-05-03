library('Meiosis')
context('meiosis')

data('exdat', package = 'Meiosis')

test_that(".meiosis_geno_", {
  opt <- list(c(1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L),
              c(1L, 1L, 1L, 1L, 1L, 0L, 1L, 0L, 0L))

  patalle <- exdat$ind$paternal[[1L]]
  matalle <- exdat$ind$maternal[[1L]]
  pos <- exdat$positions[[1L]]
  xlocations <- c(1, 18, 85, 145)

  tmp <- replicate(n = 1000L,
                   Meiosis:::.meiosis_geno_(patalle = patalle, matalle = matalle,
                                            xlocations = xlocations, pos = pos,
                                            use_names = FALSE),
                   simplify = FALSE
                   )
  tmp <- do.call(what = rbind, args = tmp)
  tmp <- unique(x = tmp)
  tmp <- split(x = tmp, f = row(tmp))
  testthat::expect_true(object = setequal(x = tmp, y = opt))
})
