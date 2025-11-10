# Adapt the number of particles

This function takes the provided \[rbi::libbi()\] and runs MCMC at a
single point (i.e., repeatedly proposing the same parameters), adapting
the number of particles distribution until the variance of the
log-likelihood crosses the value given as `target.variance` (1 by
default).

## Usage

``` r
adapt_particles(
  x,
  min = 1,
  max = 1024,
  target_variance = 1,
  quiet = FALSE,
  target.variance,
  ...
)
```

## Arguments

- x:

  a \[rbi::libbi()\] object

- min:

  minimum number of particles

- max:

  maximum number of particles

- target_variance:

  target log-likelihood variance; once this is crossed, the current
  number of particles will be used

- quiet:

  if set to TRUE, will not provide running output of particle numbers
  tested

- target.variance:

  deprecated; use `target_variance` instead

- ...:

  parameters for libbi\$run

## Value

a \[rbi::libbi()\] with the desired proposal distribution

## Examples

``` r
example_obs <- rbi::bi_read(system.file(package="rbi", "example_dataset.nc"))
example_model <- rbi::bi_model(system.file(package="rbi", "PZ.bi"))
example_bi <- rbi::libbi(model = example_model, obs = example_obs)
obs_states <- rbi::var_names(example_model, type = "obs")
max_time <- max(vapply(example_obs[obs_states], function(x) {
  max(x[["time"]])
}, 0))
if (FALSE) { # \dontrun{
  adapted <- adapt_particles(example_bi, nsamples = 128, end_time = max_time)
} # }
```
