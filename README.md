rbi.helpers
=============

[![Build Status](https://travis-ci.org/sbfnk/rbi.helpers.png?branch=master)](https://travis-ci.org/sbfnk/rbi.helpers)
  
[rbi.helpers] (https://github.com/sbfnk/RBi.helpers) is collection of helper functions to use with [rbi] (https://github.com/libbi/RBi), an R interface to [libbi] (https://github.com/libbi/LibBi), a library for Bayesian Inference.

It contains:
- `acceptance_rate`, to calculate the acceptance rate of a pMCMC run
- `sample_observations`, to sample observations from trajectories
- `adapt_mcmc`, to adapt the proposal distribution of a model according to the empirical standard deviations of accepted parameters
- `adapt_particles`, to adapt the number of particles at a given point in parameter space.
- `recreate_libbi`, to create a libbi object from past runs
- `compute_DIC`, to compute the DIC of a pMCMC run
- `plot_libbi`, to plot libbi output

Installation
=============

The easiest way to install `rbi.helpers` is to use the `devtools` package:

```r
# install.packages("devtools")
library(devtools)
install_github("sbfnk/rbi.helpers")
```
