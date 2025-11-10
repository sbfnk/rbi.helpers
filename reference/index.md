# Package index

## Adaptation functions

Functions for tuning the proposal distribution and partilce numbers.

- [`adapt_particles()`](https://sbfnk.github.io/rbi/reference/adapt_particles.md)
  : Adapt the number of particles
- [`adapt_proposal()`](https://sbfnk.github.io/rbi/reference/adapt_proposal.md)
  : Adapt the proposal distribution of MCMC using the covariance of
  samples

## Analysis functions

Functions for analysing LibBi results.

- [`acceptance_rate()`](https://sbfnk.github.io/rbi/reference/acceptance_rate.md)
  : Compute acceptance rate
- [`DIC(`*`<libbi>`*`)`](https://sbfnk.github.io/rbi/reference/DIC.md) :
  Compute Deviance Information Criterion (DIC) for a libbi model

## Date conversion

Functions for converting between libbi times and dates.

- [`numeric_to_time()`](https://sbfnk.github.io/rbi/reference/numeric_to_time.md)
  : Convert numeric times to actual times or dates
- [`time_to_numeric()`](https://sbfnk.github.io/rbi/reference/time_to_numeric.md)
  : Convert actual times or dates to numeric times
