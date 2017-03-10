library('purrr')
library('Meiosis')
## library('microbenchmark')
## bn <- microbenchmark::microbenchmark
## ls.str(as.environment("package:Meiosis"))
sim_geno_individual <- function(n_loci, alleles) {
  purrr::rerun(2L, purrr::map(n_loci, ~sample(x = alleles, size = .x, replace = TRUE)))
}


bcgv_r <- function(individual, positions, xodat_param, locus_effects,
                n_gam, se_level, min_rep = 2L, max_rep) {
  L <- vapply(X = positions, FUN = function(x) ceiling(max(x)), FUN.VALUE = numeric(1L))
  xodat_param <- Meiosis::create_xodat_param(L = L, m = 0L, p = 0.0, obligate_chiasma = FALSE)
  f <- Meiosis::meiosis_geno

  n <- 0L
  online_mean <- 0.0
  M2 <- 0.0

  vec <-numeric()
  while(n <= min_rep || (se >= se_level && n <= max_rep))
  {
    values <- numeric(n_gam)
    for (i in seq_len(n_gam)) {
      gam <- f(individual = individual, positions = positions,
                xodat_param = xodat_param)
      values[i] <- sum(vapply(X = seq_along(gam), FUN.VALUE = numeric(1L),
                              FUN = function(i) crossprod(gam[[i]], locus_effects[[i]])))
    }
    max_val <- max(values)

    vec <- c(vec,max_val)

    n <- n + 1L
    delta <- max_val - online_mean
    online_mean <- online_mean + delta / n
    delta2 <- max_val - online_mean
    M2 <- M2 + delta * delta2
    se <- sqrt(M2 / (n * (n - 1)))
  }
  return(list('bcgv' = online_mean, 'se' = se, 'n' = n))
}


n_chr <- 10L
L_ <- 300.0
m <- 0L
p <- 0.0
obligate_chiasma <- FALSE
L <- rep(x = L_, times = n_chr)
## Lstar <- purrr::map_dbl(.x = L, .f = ~Meiosis::calc_Lstar(.x, m, p, 1e-6))
xodat_param <- create_xodat_param(L = L, m = m, p = p, obligate_chiasma = obligate_chiasma)

n_loci <- rep(x = 100L, times = n_chr)
alleles_geno <- c(-1L, 1L)
positions <- purrr::map2(n_loci, L, ~sort(runif(.x, min = 0.0, max = .y)))

f1 <- create_xodat_founder(alleles = c(0L, 1L), L)
f2 <- create_xodat_founder(alleles = c(2L, 3L), L)

f1g <- sim_geno_individual(n_loci, alleles_geno)
f2g <- sim_geno_individual(n_loci, alleles_geno)

cv <- new(Meiosis::Converter, positions)
cv$insert_founder(c(0L, 1L), f1g)
cv$insert_founder(c(2L, 3L), f2g)

microbenchmark::microbenchmark(times = 1000,
x <- cross_xodat(f1, f2, xodat_param),
## cv$convert(x),
cross_geno(f1g, f2g, positions, xodat_param),
dh_geno(f1g, positions, xodat_param)
)

n_gam <- 5L
individual <- f1g
locus_effects <- purrr::map(n_loci, .f = rnorm)
min_rep <- 10L
max_rep <- 100000L
se_level <- 0.1

Meiosis::bcgv(f1g, positions, locus_effects,
              n_gam, se_level, min_rep, max_rep,
              0L, 0.0)

microbenchmark::microbenchmark(times = 1,

bcgv_r(f1g, positions, xodat_param, locus_effects,
       n_gam, se_level, min_rep, max_rep),

Meiosis::bcgv(f1g, positions, locus_effects,
              n_gam, se_level, min_rep, max_rep)
)


Meiosis::bcgv(f1g, positions, locus_effects,
              n_gam, se_level = 0.8, min_rep = 10, max_rep = 100)


cross_geno(f1g, f2g, pos, xodat_param)
cross_geno(f1g, f1g, pos, xodat_param)
dh_geno(f1g, pos, xodat_param)
self_geno(f2g, pos, xodat_param)
Meiosis::







ind <- sim_xodat_individual(rep(20L, length(L)), c(1L, 2L), L)

microbenchmark::microbenchmark(
  Meiosis::meiosis_xodat(ind, L, m, p, obligate_chiasma, Lstar),
  Meiosis::test(ind, xo_param)
)


ind <- sim_geno_individual(n_loci, alleles_geno)
bn(times = 1e2,
   Meiosis:::meiosis_geno(ind, pos, L, m, p, obligate_chiasma, Lstar),
   Meiosis:::meiosis_geno_std(ind, pos, L, m, p, obligate_chiasma, Lstar)
)

patalle <- ind[[1]][[1]]
matalle <- ind[[2]][[1]]
pos_ <- pos[[1]]
L_ <- Lstar_ <- L[[1]]

bn(times = 1e3,
Meiosis:::.meiosis_geno_std(patalle, matalle, pos_, L_, m, p, obligate_chiasma, Lstar_),
Meiosis:::.meiosis_geno(patalle, matalle, pos_, L_, m, p, obligate_chiasma, Lstar_)
)


## bn(times = 1e4,
##    meiosis_geno(ind, pos, L, m, p, obligate_chiasma, Lstar),
##    dh_geno(ind, pos, L, m, p, obligate_chiasma, Lstar),
##    cross_geno(ind, ind, pos, L, m, p, obligate_chiasma, Lstar)
## )

## l <- purrr::rerun(1000L, dh_geno(ind, pos, L, m, p, obligate_chiasma, Lstar))
## x <- do.call(rbind, at_depth(at_depth(l, 2L, unlist), 1L, ~do.call(rbind, .x)))

######################################################

sim_xodat_individual <- function(n_loci, alleles, L)
{
  purrr::rerun(2L,
    purrr::map2(n_loci, L, .f = function(.x, .y) {
      list(alleles = sample(x = alleles, size = .x, replace = TRUE),
           locations = c(sort(runif(.x - 1L, min = 0.0, max = .y)), max(.y))
           )
    })
  )
}

alleles_xodat <- c(1L, 2L)
ind <- sim_xodat_individual(rep(20L, length(L)), alleles_xodat, L)
patalle <- ind[[1]][[1]][["alleles"]]
patloc <- ind[[1]][[1]][["locations"]]
matalle <- ind[[2]][[1]][["alleles"]]
matloc <- ind[[2]][[1]][["locations"]]

microbenchmark::microbenchmark(
Meiosis::realized_ibd(ind,ind),
Meiosis::realized_f(ind)
)

ind2 <- sim_xodat_individual(n_loci, alleles_xodat, L)

attr(ind[[1]], 'ped') <- c('1' = 0.5, '3' = 0.5)
attr(ind[[2]], 'ped') <- c('1' = 0.25, '2' = 0.75)

Meiosis::meiosis_xodat(ind, L, m, p, obligate_chiasma, Lstar)

Meiosis::meiosis_xodat(ind, L, m, p, FALSE, L)
ind

attr(ind2[[1]], 'ped') <- list('b' = 1.0)
attr(ind2[[2]], 'ped') <- list('b' = 1.0)

ind1


pat_ped <- list('1' = 0.5, '3' = 0.5)
mat_ped <- list('1' = 0.25, '2' = 0.75)

alleles <- union(names(pat_ped), names(mat_ped))
ped <- list()
for (a in alleles) {
  tmp <- 0.0
  if (a %in% names(pat_ped))
    tmp <- tmp +pat_ped[[a]]

  if (a %in% names(mat_ped))
    tmp <- tmp + mat_ped[[a]]
  ped[[a]] <- tmp / 2
}

pat_ped['4'] + 8


bn(times = 1000,
  Meiosis:::.meiosis_xodat(patalle,patloc,matalle,matloc,L_,m,p,obligate_chiasma,L_),
  Meiosis:::.meiosis_xodat_std(patalle,patloc,matalle,matloc,L_,m,p,obligate_chiasma,L_)
)

microbenchmark::microbenchmark(times = 10,
                               Meiosis::realized_ibd(ind, ind2)
)

ind

f1g <- purrr::map(n_loci, ~sample(x = c(0L, 1L), size = .x, replace = TRUE))
f2g <- purrr::map(n_loci, ~sample(x = c(0L, 1L), size = .x, replace = TRUE))
f <- new(Founders)
f$insert(1, f1g)
f$insert(2, f2g)
f$size()

str(u)

microbenchmark::microbenchmark(times = 1000,
Meiosis::convert_xodat(ind, pos, f)
)



n <- 1000L
L <- 200.0
m <- 0L
p <- 0.0
obligate_chiasma <- FALSE
Lstar <- Meiosis::calc_Lstar(L, m, p, 1e-6)

f1 <- list('alleles' = 1L, 'locations' = L)
f2 <- list('alleles' = 2L, 'locations' = L)


p <- Meiosis::meiosis_xodat(f1$alleles, f1$locations, f2$alleles, f2$locations,
                            L, m, p, obligate_chiasma, Lstar)

f1g <- sample(x = c(0L, 1L), size = n, replace = TRUE)
f2g <- sample(x = c(0L, 1L), size = n, replace = TRUE)
pos <- sort(runif(n = n, min = 0, max = L))

f <- new(Founders)
f$insert(1, f1g)
f$insert(2, f2g)
f$size()

microbenchmark::microbenchmark(times = 1000,
Meiosis::convert(p$alleles, p$locations, pos, f)
)

find <- list(fgam1, fgam2)

f <- new(Founders)
u <- new(Foo, f)
u$bar(f)




geno_ind <- sim_geno_individual(n_loci, alleles_geno)
xodat_ind <- sim_xodat_individual(n_loci, alleles_xodat, L)
str(xodat_ind)


## patalle <- xodat_ind[[1]][[1]][["alleles"]]
## patloc <- xodat_ind[[1]][[1]][["locations"]]
## matalle <- xodat_ind[[1]][[2]][["alleles"]]
## matloc <- xodat_ind[[1]][[2]][["locations"]]
## xlocations <- Meiosis::crossover(L[1], m, p, obligate_chiasma, Lstar[1])

factory <- function(L, m, p, obligate_chiasma, Lstar)
{
  function(individual) {
    Meiosis::meiosis_xodat_ind(individual, L, m, p, obligate_chiasma, Lstar)
  }
}

mxiw <- factory(L, m, p, obligate_chiasma, Lstar)


microbenchmark::microbenchmark(times = 1e4L, 
Meiosis::meiosis_xodat_ind(xodat_ind, L, m, p, obligate_chiasma, Lstar),
Meiosis::meiosis_geno_ind(geno_ind, pos, L, m, p, obligate_chiasma, Lstar),
mxiw(xodat_ind)
)

Meiosis::meiosis_geno_ind(geno_ind[1], pos, L, m, p, obligate_chiasma, Lstar)

Meiosis::meiosis_xodat_ind(xodat_ind[1], L, m, p, obligate_chiasma, Lstar)

g2 <- Meiosis::meiosis_xodat_ind(xodat_ind, L, m, p, obligate_chiasma, Lstar)
g3 <- Meiosis::meiosis_xodat_ind(list(g1,g2), L, m, p, obligate_chiasma, Lstar)

g1 ms
str()
str(geno_ind)
str(gam)



n_xo <- 100L
L <- 200
m <- 0L
p <- 0.0
obligate_chiasma <- FALSE
epsilon <- 1e-6
alleles <- c(1L, 2L, 3L, 4L)

patalle <- sample(x = alleles, size = n_xo, replace = TRUE)
matalle <- sample(x = alleles, size = n_xo, replace = TRUE)
patloc <- c(sort(runif(n = n_xo - 1L, min = 0, max = L)), L)
matloc <- c(sort(runif(n = n_xo - 1L, min = 0, max = L)), L)
xo <- Meiosis::crossover(L, m, p, obligate_chiasma, L)
Meiosis::meiosis_xodat(patalle, patloc, matalle, matloc, xo)


paternal <- purrr::rerun(n_chr, list(patalle, patloc))
maternal <- purrr::rerun(n_chr, list(matalle, patloc))
individual <- list(paternal, maternal)

pos <- purrr::rerun(n_chr, sort(runif(n = n_loci, min = 0, max = L)))

gam1 <- meiosis_individual(individual, Lvec, m = 0L, p = 1.0,  FALSE, Lvec)
gam2 <- meiosis_individual(individual, Lvec, m = 0L, p = 1.0,  FALSE, Lvec)
ind <- list(gam1, gam2)


microbenchmark(times = 1e4L, 
meiosis_individual(ind, Lvec, m = 0L, p = 1.0,  FALSE, Lvec),
test(individual =ind, L = Lvec, m = 0L, p = 1.0,  FALSE, Lvec)
)


test <- function(individual, L, m, p, obligate_chiasma, Lstar)
{
  n_chr <- length(L)
  gamete <- vector('list', n_chr)
  paternal <- individual[[1]]
  maternal <- individual[[2]]
  for (i in seq_len(n_chr)) {
    pat <- paternal[[i]]
    mat <- maternal[[i]]
    gamete[[i]] <- Meiosis::meiosis_xodat_R(pat[[1]], pat[[2]], mat[[1]], mat[[2]],
                      Meiosis::crossover(L[[i]], m, p, obligate_chiasma, Lstar[[i]]))
  }
  gamete
}



patalle_ <- sample(x = alleles, size = n_loci, replace = TRUE)
matalle_ <- sample(x = alleles, size = n_loci, replace = TRUE)

xo <- Meiosis::crossover(L = L[1], m = 0L, p = 1.0, obligate_chiasma = FALSE, Lstar = L[1])
Meiosis::meiosis_xodat_R(patalle, patloc, matalle, matloc, xo)


founder <- purrr::rerun(.n = length(alleles), sample(x = c(0L, 1L), size = n_loci, replace = TRUE))
names(founder) <- alleles

pos <- sort(runif(n = n_loci, min = 0, max = L))

microbenchmark(times = 1e4L, 
a = {
  xo <- Meiosis::crossover(L = L, m = 0L, p = 1.0, obligate_chiasma = FALSE, Lstar = L)
  xodat <- Meiosis::meiosis_xodat(patalle, patloc, matalle, matloc, xo)
  Meiosis::convert(xodat$alleles, xodat$locations, pos, founder)
},
b = {
  xo <- Meiosis::crossover(L = L, m = 0L, p = 1.0, obligate_chiasma = FALSE, Lstar = L)
  Meiosis::meiosis_geno(patalle_, matalle_, xo, pos)
}
)

patalle_ <- sample(x = alleles, size = n_loci, replace = TRUE)
matalle_ <- sample(x = alleles, size = n_loci, replace = TRUE)

microbenchmark(times = 1e4L, 
               xo <- Meiosis::crossover(L = L, m = 0L, p = 1.0, obligate_chiasma = FALSE, Lstar = L),
               Meiosis::meiosis_xodat(patalle, patloc, matalle, matloc, xo),
               Meiosis::meiosis_geno(patalle_, matalle_, xo, pos)
)

