Helper functions for the rbi package
================

<!-- badges: start -->

![GitHub R package
version](https://img.shields.io/github/r-package/v/sbfnk/rbi.helpers)
[![R-CMD-check](https://github.com/sbfnk/rbi/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sbfnk/rbi/actions/workflows/R-CMD-check.yaml)
[![codecov](https://app.codecov.io/github/sbfnk/rbi.helpers)](https://codecov.io/github/sbfnk/rbi.helpers)
![GitHub
contributors](https://img.shields.io/github/contributors/sbfnk/rbi.helpers)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN
status](https://www.r-pkg.org/badges/version/rbi.helpers)](https://CRAN.R-project.org/package=rbi.helpers)
<!-- badges: end -->

[rbi.helpers](https://github.com/sbfnk/rbi.helpers) is collection of
helper functions to use with [rbi](https://github.com/sbfnk/rbi), an R
interface to [LibBi](https://github.com/lawmurray/LibBi), a library for
Bayesian Inference.

It contains:

- `adapt_proposal`, to adapt the proposal distribution of a model
  according to the empirical standard deviations of accepted parameters
- `adapt_particles`, to adapt the number of particles at a given point
  in parameter space.
- `DIC`, to compute the DIC of a pMCMC run
- `acceptance_rate`, to calculate the acceptance rate of a pMCMC run

# Installation

The easiest way to install `rbi.helpers` is to use the `remotes`
package:

``` r
# install.packages("remotes")
remotes::install_github("sbfnk/rbi.helpers")
```
