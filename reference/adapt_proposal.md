# Adapt the proposal distribution of MCMC using the covariance of samples

This function takes the provided \[rbi::libbi()\] object and runs MCMC,
adapting the proposal distribution until the desired acceptance rate is
achieved. If a scale is given, it will be used to adapt the proposal at
each iteration

## Usage

``` r
adapt_proposal(
  x,
  min = 0,
  max = 1,
  scale = 2,
  max_iter = 10,
  adapt = c("size", "shape", "both"),
  size = FALSE,
  correlations = TRUE,
  truncate = TRUE,
  quiet = FALSE,
  ...
)
```

## Arguments

- x:

  `link{libbi}` object

- min:

  minimum acceptance rate

- max:

  maximum acceptance rate

- scale:

  scale multiplier/divider for the proposal. If \>1 this will be
  inverted.

- max_iter:

  maximum of iterations (default: 10)

- adapt:

  what to adapt; if "size" (default), the width of independent proposals
  will be adapted; if "shape", proposals will be dependent (following a
  multivariate normal) taking into account empirical correlations; if
  "both", the size will be adapted before the shape

- size:

  (deprecated, use `{adapt}` instead) if TRUE (default: FALSE), the size
  of the (diagonal multivariate normal) proposal distribution will be
  adapted

- correlations:

  (deprecated, use `{adapt}` instead) if TRUE (default: FALSE), the
  shape of the (diagonal multivariate normal) proposal distribution will
  be adapted according to the empirical covariance

- truncate:

  if TRUE, the proposal distributions will be truncated according to the
  support of the prior distributions

- quiet:

  if set to TRUE, will not provide running output of particle numbers
  tested

- ...:

  parameters for \[rbi::sample()\]

## Value

a \[rbi::libbi()\] with the desired proposal distribution

## Examples

``` r
example_obs <- rbi::bi_read(system.file(package="rbi", "example_dataset.nc"))
example_model <- rbi::bi_model(system.file(package="rbi", "PZ.bi"))
example_bi <- rbi::libbi(model = example_model, obs = example_obs)
obs_states <- rbi::var_names(example_model, type="obs")
max_time <- max(vapply(example_obs[obs_states], function(x) {
  max(x[["time"]])
}, 0))
# adapt to acceptance rate between 0.1 and 0.5
if (FALSE) { # \dontrun{
  adapted <- adapt_proposal(example_bi,
    nsamples = 100, end_time = max_time,
    min = 0.1, max = 0.5, nparticles = 256, correlations = TRUE
  )
} # }
```
