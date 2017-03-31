Meiosis: Simulation of meiosis in plant breeding research
======
[![Travis-CI Build Status](https://travis-ci.org/DominikMueller64/Meiosis.svg?branch=master)](https://travis-ci.org/DominikMueller64/Meiosis)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/Meiosis)](https://cran.r-project.org/package=Meiosis)
[![Coverage Status](https://img.shields.io/codecov/c/github/DominikMueller64/Meiosis/master.svg)](https://codecov.io/github/DominikMueller64/Meiosis?branch=master)

**Meiosis** is a lean [R](http://www.r-project.org) package for the simulation
of meiosis events in diploid (or allo-polyploid) plant species for genetic
research in plant breeding.

---

### Installation

You can install bcgvr from its [GitHub repository](http://github.com/DominikMueller64/Meiosis).
You first need to install the [devtools](https://github.com/hadley/devtools) package.

```r
install.packages("devtools")
```

Then install Meiosis with 

```r
devtools::install_github("DominikMueller64/Meiosis", build_vignettes = TRUE)
```

Windows users need to make sure that [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
is installed.

---

### Example

A simple example for simulating Meiosis based on genotypic data is give here:

```r
## Simulate some data
L <- c(332, 221) ## length of chromosomes
n_loci <- c(20L, 43L) ## number of loci on chromosomes
## Simulate genetic positions of loci
positions <- lapply(seq_along(n_loci), function(i) sort(runif(n_loci[i], 0, L[i])))
## Simulate genotypic data for an individual
ind <- replicate(2L, lapply(n_loci, function(n) sample(c(0L, 1L), n, TRUE)), simplify = FALSE)
## Construct a parameter list necessary for simulating Meiosis (needs to be done once)
xoparam <- Meiosis::create_xoparam(L)

Meiosis::meiosis_geno(ind, positions, xoparam) ## Simulate a new gamete
Meiosis::cross_geno(ind, ind, positions, xoparam) ## Simulate a new individual
Meiosis::self_geno(ind, positions, xoparam) ## Simulate a new selfing
Meiosis::dh_geno(ind, positions, xoparam) ## Simulate a new doubled haploid
```

---

### Vignette

A vignette describing the functionality of the package and the data structures
is available from within R. 

```r
vignette('Meiosis', package = 'Meiosis')
```

---

### Acknowledgements

Parts of the core functionality and documentation of **Meiosis** was inspired and
adapted, respectively, from the package
[simcross](https://github.com/kbroman/simcross) of [Karl
Broman](http://kbroman.org/).

---

### Author

Dominik Mueller

---

### License

This package is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License, version 3, as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License for more details.

A copy of the GNU General Public License, version 3, is available at
<http://www.r-project.org/Licenses/GPL-3>

