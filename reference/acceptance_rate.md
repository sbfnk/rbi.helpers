# Compute acceptance rate

This function takes the provided \[rbi::libbi()\] object which has been
run, or a bi file, and returns a the acceptance rate

## Usage

``` r
acceptance_rate(...)
```

## Arguments

- ...:

  parameters to \[rbi::get_traces()\] (especially 'x')

## Value

acceptance rate

## Author

Sebastian Funk

## Examples

``` r
example_run <- rbi::bi_read(
  system.file(package = "rbi.helpers", "example_run.nc")
)
example_model_file <- system.file(package = "rbi", "PZ.bi")
example_bi <- rbi::attach_data(
  rbi::libbi(example_model_file), "output", example_run
)
acceptance_rate(example_bi)
#> [1] 0.4040404
```
