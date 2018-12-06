Helper functions for the RBi package
=============

[![Build Status](https://travis-ci.org/sbfnk/RBi.helpers.png?branch=master)](https://travis-ci.org/sbfnk/RBi.helpers) [![codecov](https://codecov.io/github/sbfnk/RBi/branch/master/graphs/badge.svg)](https://codecov.io/github/sbfnk/RBi) 
  
[RBi.helpers](https://github.com/sbfnk/RBi.helpers) is collection of helper functions to use with [RBi](https://github.com/libbi/RBi), an R interface to [LibBi](https://github.com/libbi/LibBi), a library for Bayesian Inference.

It contains:
- `adapt_proposal`, to adapt the proposal distribution of a model according to the empirical standard deviations of accepted parameters
- `adapt_particles`, to adapt the number of particles at a given point in parameter space.
- `DIC`, to compute the DIC of a pMCMC run
- `acceptance_rate`, to calculate the acceptance rate of a pMCMC run

Installation
=============

The easiest way to install `rbi.helpers` is to use the `devtools` package:

```r
# install.packages("devtools")
library(devtools)
install_github("sbfnk/rbi.helpers")
```
