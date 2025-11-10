# Compute Deviance Information Criterion (DIC) for a libbi model

Computes the DIC of a libbi object containing Monte-Carlo samples. The
effective number of parameters is calculated following Gelman et al.,
Bayesian Data Analysis: Second Edition, 2004, p. 182.

## Usage

``` r
# S3 method for class 'libbi'
DIC(x, bootstrap = 0, ...)
```

## Arguments

- x:

  a `libbi` object

- bootstrap:

  number of bootstrap samples to take, 0 to just take data

- ...:

  any parameters to be passed to \`bi_read\` (e.g., \`burn\`)

## Value

DIC

## Author

Sebastian Funk

## Examples

``` r
example_run <- rbi::bi_read(
  system.file(package = "rbi", "example_output.nc")
)
example_model_file <- system.file(package = "rbi", "PZ.bi")
example_bi <- rbi::attach_data(
  rbi::libbi(example_model_file), "output", example_run
)
DIC(example_bi)
#> [1] 173.834
```
