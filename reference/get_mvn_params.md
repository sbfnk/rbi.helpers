# Construct a covariance matrix

This function takes the provided \[rbi::libbi()\] which has been run and
returns the square root of the covariance matrix, which can be used for
proposal distributions

## Usage

``` r
get_mvn_params(x, scale = 1, correlations = TRUE, fix)
```

## Arguments

- x:

  a \[rbi::libbi()\] which has been run

- scale:

  a factor by which to scale all elements of the covariance matrix

- correlations:

  logical; if TRUE, correlations are taken into account when
  constructing the parameters

- fix:

  if this is set, all elements of the covariance matrix will be set to
  it

## Value

the updated bi model
