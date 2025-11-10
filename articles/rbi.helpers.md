# Introduction to rbi.helpers

[rbi.helpers](https://github.com/sbfnk/rbi.helpers) is collection of
helper functions to use with [rbi](https://github.com/sbfnk/rbi), an R
interface to [LibBi](https://github.com/lawmurray/LibBi), a library for
Bayesian Inference.

This vignette builds on the [rbi
vignette](https://sbfnk.github.io/rbi/articles/introduction.html),
applying the higher-level functions contained in **rbi.helpers** to the
same model introduced there. For the lower-level functions to run
**LibBi** and read the results, please refer to the documentation and
vignette that comes with **rbi**.

## Installation

The easiest way to install the latest stable version of **rbi.helpers**
is via CRAN. The package name is `rbi.helpers` (all lowercase):

``` r
install.packages("rbi.helpers")
```

Alternatively, the current development version can be installed using
the `remotes` package

``` r
remotes::install_github("sbfnk/rbi.helpers")
```

Most functions in the **rbi.helpers** package require a working
installation of **LibBi**. Please see the [rbi
vignette](https://sbfnk.github.io/rbi/articles/introduction.html) for
how to get one via homebrew or linuxbrew.

## Loading the package

Use

``` r
library("rbi")
library("rbi.helpers")
```

to load the package.

## Loading the model and generating a synthetic dataset

These steps are reproduced from the [rbi
vignette](https://sbfnk.github.io/rbi/articles/introduction.html), where
there is more information on the individual steps

``` r
model_file <- system.file(package = "rbi", "SIR.bi") # file included in package
sir_model <- bi_model(model_file) # load model
set.seed(1001912)
sir_data <- bi_generate_dataset(sir_model, end_time = 16 * 7, noutputs = 16)
#> Warning in bi_generate_dataset(sir_model, end_time = 16 * 7, noutputs = 16):
#> `bi_generate_dataset` has been deprecated and will be removed from future
#> versions of `rbi`. Use the `generate_dataset` function instead.
```

## Adapt the number of particles

In the [rbi
vignette](https://sbfnk.github.io/rbi/articles/introduction.html), a
[stochastic SIR
model](https://raw.githubusercontent.com/sbfnk/rbi/master/inst/SIR.bi)
is fitted to simulated data from the same model using particle
Markov-chain Monte Carlo with 32 particles. Given a model and data, how
do we know how many particles we need? This question does not have a
simple answer, as the “optimal” number of particles may depend on the
state of the Markov chain. A possible rule-of-thumb is to choose the
number of particles such that the variance of the log-likelihood near
the mode is approximately one. This suggests a
[strategy](https://darrenjw.wordpress.com/2014/06/08/tuning-particle-mcmc-algorithms/)
by which first and approximate location of the mode or mean of the
posterior distribution is obtained in a trial run, before the numer of
particles is adjusted by monitoring the variance of the log-likelihood
while keeping the parameters fixed. **rbi.helpers** implements the
second part of this strategy (adjusting the number of particles at a
given location in parameter space) with the `adapt_particles` method.
For the first part (finding the mode), a crude method is to take a fixed
number of samples from the prior distribution and choose the one that
maximises the posterior. In **rbi**, this can be achieved with

``` r
bi_prior <- sample(
  proposal = "prior", sir_model, nsamples = 1000, end_time = 16 * 7,
  nparticles = 4, obs = sir_data, seed = 1234
)
```

This runs particle MCMC with the prior distribution as proposal
distribution (because `proposal = "prior"`), in this case with 4
particles. In other words, when sampling from the posterior the
proposals will be drawn independently from the prior distribution. Note
that we set a seed to make the results reproducible. It is worth trying
the commands with a different seed and seeing the difference to the
results obtained below. The location in parameters of the sampler at the
end of the 1000 iterations will give an approximation of the mode of the
posterior distribution. This can then be used to adjust the number of
particles using

``` r
adapted <- adapt_particles(bi_prior)
#> Mon Nov 10 09:12:29 2025 Adapting the number of particles
#> Mon Nov 10 09:12:45 2025 4 particles, loglikelihod variance: 9.5863402696008
#> Mon Nov 10 09:12:49 2025 8 particles, loglikelihod variance: 3.77192932450897
#> Mon Nov 10 09:12:53 2025 16 particles, loglikelihod variance: 2.24273655866084
#> Mon Nov 10 09:12:57 2025 32 particles, loglikelihod variance: 1.30179246083847
#> Mon Nov 10 09:13:05 2025 64 particles, loglikelihod variance: 0.651044197683024
#> Mon Nov 10 09:13:05 2025 Choosing 64 particles.
```

This will take the last sample of the output file contained in the
`libbi` object `bi_prior`, and use it to adjust the number of particles
by starting with 1 particle (or a given `min`) and doubling it until the
variance of the loglikelihood crosses 1. The number of particles is then
saved in the `adapted` object:

``` r
adapted$options$nparticles
#> [1] 64
```

## Adapt the proposal distribution

Having adjusted the number of particles, the second important
information to give the posterior sampler is the proposal distribution.
This can, again, be obtained using a sequence of trial runs, whereby the
proposal distribution is sequentially adjusted from previous samples to
be proportional to the empirical covariance of the posterior samples.
The way this is implemented in the `adapt_proposal` function in
**rbi.helpers** is that the proposal distribution is adjusted to come
from a multivariate normal taking into account the covariance of samples
obtained so far, until the acceptance rate lies between the given
minimum and maximumad. For example, to adjust the proposal distribution
for an acceptance rate between 0.05 and 0.4, we can run:

``` r
adapted <- adapt_proposal(adapted, min = 0.05, max = 0.4)
#> Mon Nov 10 09:13:05 2025 Adapting the proposal distribution
#> Mon Nov 10 09:13:05 2025 Initial trial run
#> Mon Nov 10 09:13:26 2025 Acceptance rate: 0.283283283283283
```

The covariance matrices for parameters and initial conditions are stored
in the input file contained in the `libbi` object `adapted`

``` r
bi_read(adapted, file = "input")
#> $`__proposal_parameter_cov`
#>   __dim_parameter_cov.1 __dim_parameter_cov.2      value
#> 1                 p_rep                 p_rep 0.03184091
#> 2                 p_rep                  p_R0 0.00000000
#> 3                  p_R0                 p_rep 0.00000000
#> 4                  p_R0                  p_R0 0.25829530
```

## Compute DIC

To compute the [Deviance Information
Criterion](https://en.m.wikipedia.org/wiki/Deviance_information_criterion)
(DIC), use `DIC`:

``` r
posterior <- sample(adapted)
DIC(posterior)
#> [1] 72.00597
```

This samples from the posterior distribution using the adapted number of
particles and proposal distribution and then calculates the DIC.

## Convert between LibBi times and actual times or dates

LibBi uses real-valued times. To convert these to time or date objects
for use in R, use the `numeric_to_time` function:

``` r
res <- numeric_to_time(posterior, unit = "day", origin = as.Date("2018-04-01"))
head(res$Z)
#>   np       time      value
#> 1  0 2018-04-01   1.000000
#> 2  0 2018-04-08   1.452416
#> 3  0 2018-04-15   6.895681
#> 4  0 2018-04-22  17.104379
#> 5  0 2018-04-29 134.397199
#> 6  0 2018-05-06 365.274828
```

The function `time_to_numeric` does the converse, converting R times or
dates into numeric values for use by LibBi:

``` r
orig <- time_to_numeric(res, unit = "day", origin = as.Date("2018-04-01"))
head(orig$Z)
#>   np time      value
#> 1  0    0   1.000000
#> 2  0    7   1.452416
#> 3  0   14   6.895681
#> 4  0   21  17.104379
#> 5  0   28 134.397199
#> 6  0   35 365.274828
```

## Create inference chains

With the pipe operator available since R version 4.1, it is possible to
construct inference chains more elegantly. For example, to get posterior
samples from adapted Metropolis-Hastings including sampled observations,
we could have written

``` r
posterior <- sample(
  proposal = "prior", sir_model, nsamples = 1000,
  end_time = 16 * 7, nparticles = 4, obs = sir_data, seed = 1234
) |>
  adapt_particles() |>
  adapt_proposal(min = 0.05, max = 0.4) |>
  sample(nsamples = 5000) |>
  sample_obs()
#> Mon Nov 10 09:13:47 2025 Adapting the proposal distribution
#> Mon Nov 10 09:14:15 2025 Adapting the number of particles
#> Mon Nov 10 09:14:32 2025 4 particles, loglikelihod variance: 4.98864430510319
#> Mon Nov 10 09:14:35 2025 8 particles, loglikelihod variance: 2.50334348358283
#> Mon Nov 10 09:14:41 2025 16 particles, loglikelihod variance: 2.278740962349
#> Mon Nov 10 09:14:46 2025 32 particles, loglikelihod variance: 2.06221670510984
#> Mon Nov 10 09:14:53 2025 64 particles, loglikelihod variance: 1.02667431640938
#> Mon Nov 10 09:15:05 2025 128 particles, loglikelihod variance: 0.399685704417314
#> Mon Nov 10 09:15:05 2025 Choosing 128 particles.
#> Mon Nov 10 09:15:06 2025 Initial trial run
#> Mon Nov 10 09:15:32 2025 Acceptance rate: 0.388388388388388
```
