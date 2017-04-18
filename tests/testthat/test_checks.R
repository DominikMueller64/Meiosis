library('Meiosis')
context("checks")

data('exdat', package = 'Meiosis')

x <- exdat$positions
test_that("check_positions ok", {
  testthat::expect_null(check_positions(x))
  y <- x
  y[[1L]] <- rev(y[[1L]])
  testthat::expect_error(check_positions(y))
  y <- x
  y[[1L]][[1L]] <- -1
  testthat::expect_error(check_positions(y))
  y <- x
  y[[1L]][[2L]] <- y[[1L]][[1L]]
  testthat::expect_error(check_positions(y))
})

x <- exdat$founder
test_that("check_xo_individual ok", {
  testthat::expect_null(check_xo_individual(x))
  testthat::expect_error(check_xo_individual(x[[1L]]))
  y <- x
  y[[1L]][[1L]] <- NULL
  testthat::expect_error(check_xo_individual(y))
  y <- x
  y[[1L]][[1L]][[1L]] <- as.integer(rnorm(100L))
  testthat::expect_error(check_xo_individual(y))
})

x <- exdat$ind
test_that("check_geno_individual ok", {
  testthat::expect_null(check_geno_individual(x))
  testthat::expect_error(check_geno_individual(x[[1L]]))
  x[[1L]][[1L]] <- NULL
  testthat::expect_error(check_geno_individual(x))
})

