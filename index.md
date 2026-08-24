# Installation

Helper functions for the rbi package [![rbi.helpers
website](reference/figures/logo.png)](https://sbfnk.github.io/rbi.helpers/)
================

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

The easiest way to install `rbi.helpers` is to use the `remotes`
package:

``` r

# install.packages("remotes")
remotes::install_github("sbfnk/rbi.helpers")
```
