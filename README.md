rbi.helpers
=============

[![Build Status](https://travis-ci.org/sbfnk/RBi.helpers.png?branch=master)](https://travis-ci.org/sbfnk/RBi.helpers)
  
[rbi.helpers] (https://github.com/sbfnk/RBi.helpers) is collection of helper functions to use with [rbi] (https://github.com/libbi/RBi), an R interface to [libbi] (https://github.com/libbi/LibBi), a library for Bayesian Inference.

It contains:
- `adapt_mcmc`, to adapt the proposal distribution of a model according to the empirical standard deviations of accepted parameters
- `adapt_particles`, to adapt the number of particles at a given point in parameter space.
- `plot_libbi`, to plot libbi output
- `compute_DIC`, to compute the DIC of a pMCMC run
- `acceptance_rate`, to calculate the acceptance rate of a pMCMC run
- `recreate_libbi`, to create a libbi object from past runs

Installation
=============

The easiest way to install `rbi.helpers` is to use the `devtools` package:

```r
# install.packages("devtools")
library(devtools)
install_github("sbfnk/rbi.helpers")
```
