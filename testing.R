set.seed(123L)
## The C++ routines use an independent random number generator. For seeding it, do e.g.
Meiosis::seed_rng(seed = 123L)

## Example: Create crossover-parameters
n_chr <- 3L  ## number of chromosomes
L <- runif(n = n_chr, min = 100, max = 300)  ## sample length of chromosomes in cM
xoparam <- create_xoparam(L)  ## no interference, no obligate chiasma
str(xoparam)


## Genotypic data: number of loci per chromosome and positions.
n_loci <- round(runif(n = n_chr, min = 5L, max = 10L))  ## sample number of loci per chromosome
## sample positions of loci on the chromosome
positions <- lapply(seq_len(n_chr), function(i) sort(runif(n_loci[i], min = 0, max = L[i])))

## Example 1: Simulate meiosis with genotypes.
ind <- replicate(2L, lapply(n_loci, function(n) sample(c(0L, 1L), n, replace = TRUE)),
                 simplify = FALSE) ## simulate some genotypic data
str(ind)

p_geno <- Meiosis::cross_geno(father = ind, mother = ind, positions = positions,
                              xoparam = xoparam)
str(p_geno)


## Example 2: Simulate meiosis with segmental representation.
f_alleles <- c(21L, 65L)
f <- Meiosis::create_xo_founder(alleles = f_alleles, L = L)

p_xo <- Meiosis::cross_xo(father = f, mother = f, xoparam = xoparam)
str(p_xo)

## Create a converter for converting segmental data to genotypic data.
conv <- new(Meiosis::Converter, positions)
conv$insert_founder(f_alleles, ind)
conv$convert(f)
conv$convert(f)

conv
?Meiosis::Converter
## Example 3: Derive n inbred lines from a bi-parental cross.
n_self <- 10L  ## number of generations of selfing
n <- 30L ## number of progeny

## Genotypic representation

## Second individual as parent.
ind2 <- replicate(2L, lapply(n_loci, function(n) sample(c(0L, 1L), n, replace = TRUE)),
                  simplify = FALSE) ## simulate some genotypic data

pop <- replicate(n, Meiosis::cross_geno(ind, ind2, positions, xoparam), simplify = FALSE)
for (i in seq_len(n))
  for (j in seq_len(n))
    pop[[i]] <- Meiosis::cross_geno(pop[[i]], pop[[i]], positions, xoparam)

## Segmental representation
f2 <- create_xo_founder(alleles = c(55L, 77L), L = L)
pop_xo <- replicate(n, Meiosis::cross_xo(f, f2, xoparam), simplify = FALSE)
for (i in seq_len(n))
  for (j in seq_len(n))
    pop_xo[[i]] <- Meiosis::cross_xo(pop_xo[[i]], pop_xo[[i]], xoparam)

## conv$convert(pop[[1]]) ## error, because genotypic data of second founder not present
conv$insert_founder(c(55L, 77L), ind2)  ## insert second founder first
pop_geno <- lapply(pop_xo, conv$convert) ## convert whole population

## Further examples

## Simulate a doubled haploid individual.
Meiosis::dh_geno(ind, positions, xoparam)
conv$convert(Meiosis::dh_xo(f, xoparam))

## Calculate realized coefficients of co-ancestry.
Meiosis::realized_coancestry(f)
Meiosis::realized_coancestry(p_xo) ## selfing progeny, expected coefficient of coancestry is 0.75.
Meiosis::realized_coancestry(pop_xo[[1L]], pop_xo[[2L]]) ## realized CoC of two full-sibs.

exdat <- list(xoparam = xoparam,
              positions = positions,
              ind = ind,
              ind2 = ind2,
              founder = f,
              founder2 = f2)
devtools::use_data(exdat)

library('Meiosis')

library('simcross')

x <- rnorm(10)
x <- c(1, 2, 3, 4, 4)
Meiosis::is_sorted(x, T, T)
Meiosis::check_positions(x)

positions <- list(1:10, 1:40)
conv <- new("Converter", positions)
?conv$convert
str(conv)
conv$.CppClassDef@docstring

methods::promptClass('Converter')
conv$insert_founder(c(5, 10), list(list(1:10, 1:39), list(1:10, 1:40)))

conv$size()
conv$test("bla")

conv$getRefClass()

conv$more

conv$.CppObject$size()

cv <- Meiosis::Converter$new()
cv$more
cv$doner
cv$size()
  cv
new(Meiosis::Converter)
Meiosis::Converter()
set.seed(123L)

Meiosis::seed_rng(55)
microbenchmark::microbenchmark(times = 1e4L,
     crossover(100, 100, 0.9, F, 100),
     sim_crossovers(100, 100, 0.9, F, 100)
                               )
u <- Meiosis::Converter()
u$initialize()
u$size()
u$insert_founder()

Meiosis::Converter
Meiosis::.__C__Rcpp_Converter
wtf <- Meiosis::genConverter()
wtf$doner
wtf$initialize()
wtf
new(wtf, positions)

new(Meiosis::Converter)
new(Meiosis::.__C__Converter)
new(Meiosis::.__C__Rcpp_Converter)







n_chr <- 10L
L_ <- 300.0
m <- 0L
p <- 0.0
obligate_chiasma <- FALSE
L <- rep(x = L_, times = n_chr)
## Lstar <- purrr::map_dbl(.x = L, .f = ~Meiosis::calc_Lstar(.x, m, p, 1e-6))

n_loci <- rep(x = 10000L, times = n_chr)


## n_loci <- 10000L
## L <- 500.0
## n_iter <- 1e2L

## pred <- vector('list', n_iter)
## for (i in seq_len(n_iter)) {
##   pred[[i]]$patalle <- sample(c(0L, 1L, 2L), size = n_loci, replace = TRUE)
##   pred[[i]]$matalle <- sample(c(5L, 8L, 9L), size = n_loci, replace = TRUE)
##   pred[[i]]$patloc <- c(sort(runif(n_loci - 1, max = L)), L)
##   pred[[i]]$matloc <- c(sort(runif(n_loci - 1, max = L)), L)
## }

## i <- pred[[1]]
## meiosis_xodat_test_R(i$patalle, i$patloc, i$matalle, i$matloc,L,0L,0.0,F, L)

## microbenchmark::microbenchmark(times = 100L,
## a = {for(i in pred) meiosis_xodat_test_R(i$patalle, i$patloc, i$matalle, i$matloc,L,0L,0.0,F,L)},
## b = {for(i in pred) meiosis_xodat_test_R2(i$patalle, i$patloc, i$matalle, i$matloc,L,0L,0.0,F, L)}
## )


## library('microbenchmark')
## bn <- microbenchmark::microbenchmark
## ls.str(as.environment("package:Meiosis"))
sim_geno_individual <- function(n_loci, alleles) {
  purrr::rerun(2L, purrr::map(n_loci, ~sample(x = alleles, size = .x, replace = TRUE)))
}


library("Meiosis")

n_chr <- 10L
L_ <- 300.0
m <- 0L
p <- 0.0
obligate_chiasma <- FALSE
L <- rep(x = L_, times = n_chr)
## Lstar <- purrr::map_dbl(.x = L, .f = ~Meiosis::calc_Lstar(.x, m, p, 1e-6))
xodat_param <- create_xodat_param(L = L, m = m, p = p, obligate_chiasma = obligate_chiasma)

n_loci <- rep(x = 10000L, times = n_chr)
alleles_geno <- c(-1L, 1L)
positions <- purrr::map2(n_loci, L, ~sort(runif(.x, min = 0.0, max = .y)))

n_f <- 10L
f <- lapply(seq_len(n_f), function(i)create_xodat_founder(c(2L*i-1L, 2L*i), L))
## f1 <- create_xodat_founder(alleles = c(0L, 1L), L)
## f2 <- create_xodat_founder(alleles = c(2L, 3L), L)
fg <- replicate(n_f, sim_geno_individual(n_loci, alleles_geno), FALSE)
## f1 <- create_xodat_founder(alleles = c(0L, 1L), L)

## f1g <- sim_geno_individual(n_loci, alleles_geno)
## f2g <- sim_geno_individual(n_loci, alleles_geno)

cv <- new(Meiosis::Converter, positions)
for (i in seq_len(n_f)) cv$insert_founder(c(2L*i-1L, 2L*i), fg[[i]])

## cv$insert_founder(c(0L, 1L), f1g)
## cv$insert_founder(c(2L, 3L), f2g)

cv2 <- new(Meiosis::Converter2, positions)
for (i in seq_len(n_f)) cv2$insert_founder(c(2L*i-1L, 2L*i), fg[[i]])
## cv2$insert_founder(c(0L, 1L), f1g)
## cv2$insert_founder(c(2L, 3L), f2g)

cv3 <- new(Meiosis::Converter3, positions)
for (i in seq_len(n_f)) cv3$insert_founder(c(2L*i-1L, 2L*i), fg[[i]])

x <- cross_xodat(f[[1]], f[[2]], xodat_param)

x[[1]][[1]][[1]] <- c(1,2,1,2,1)
microbenchmark::microbenchmark(times = 1000,
                               cv$convert(x, TRUE),
                               cv$convert(x, FALSE),
                               )
str(u, m =2)

microbenchmark::microbenchmark(times = 1000,
cv$convert(x),
cv2$convert(x),
cv3$convert(x),
)

all(cv3$convert(x)[[1]][[10]] == cv2$convert(x)[[1]][[10]])

u <- cv2$convert(x)
str(u, m =2)



microbenchmark::microbenchmark(times = 100,
for (i in 1:50) y <- cross_geno(f1g, f2g, positions, xodat_param),
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

