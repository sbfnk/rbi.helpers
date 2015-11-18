RBi.helpers
=============

[RBi.helpers] (https://github.com/sbfnk/RBi.helpers) is collection of helper functions to use with [RBi] (https://github.com/libbi/RBi), an R interface to [libbi] (https://github.com/libbi/LibBi), a library for Bayesian Inference.

It contains:
- `get_traces` an interface to [coda] (https://cran.r-project.org/web/packages/coda/index.html)
- `acceptance_rate`, to calculate the acceptance rate of a pMCMC run
- `sample_observations`, to sample observations from trajectories
- `adapt_mcmc`, to adapt the proposal distribution of a model according to the empirical standard deviations of accepted parameters
- `adapt_particles`, to adapt the number of particles at a given point in parameter space.
- `recreate_libbi`, to create a libbi object that from past runs
- `plot_libbi`, to plot libbi output

Installation
=============

The easiest way to install `RBi.helpers` is to use the `devtools` package:

```r
# install.packages("devtools")
library(devtools)
install_github("sbfnk/RBi.helpers")
```
