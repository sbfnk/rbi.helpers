# Construct a proposal from run results

This function takes the provided \[rbi::bi_model()\] and adds a generic
Gaussian proposal distribution.

## Usage

``` r
update_proposal(model, truncate = TRUE, blocks = c("parameter", "initial"))
```

## Arguments

- model:

  a \[rbi::bi_model()\] object

- truncate:

  truncate the multivariate normal proposals according to the used
  priors, e.g. truncating a parameter with beta prior at 0 and 1

- blocks:

  blocks to use (out of "parameter" and "initial")

## Value

the updated bi model
