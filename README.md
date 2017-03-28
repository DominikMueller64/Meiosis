

Meiosis
======
[![Travis-CI Build Status](https://travis-ci.org/DominikMueller64/Meiosis.svg?branch=master)](https://travis-ci.org/DominikMueller64/Meiosis)

Meiosis is an [R](http://www.r-project.org) package for the simulation of meiosis events.
The package serves a dual purpose. 


computation/estimation of
best conditional gametic values. The package contains a fast C++ routine and an
[Rcpp](http://www.rcpp.org/)-wrapper for exposing it to R.

[//]: # (TODO: Add reference to publication.)

---

### Installation

You can install bcgvr from its [GitHub repository](http://github.com/DominikMueller64/bcgvr).
You first need to install the [devtools](https://github.com/hadley/devtools) package.

```r
install.packages("devtools")
```

Then install bcgvr with 

```r
devtools::install_github("DominikMueller64/bcgvr")
```

Windows users need to make sure that [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
is installed.

---

### Example

Please refer to the documentation of the function `bcgv` for an example of usage.


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


