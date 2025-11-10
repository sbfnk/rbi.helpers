# Convert numeric times to actual times or dates

This function converts from numeric times (i.e., 0, 1, 2, ...) to actual
times or dates

## Usage

``` r
numeric_to_time(x, origin, unit, ...)
```

## Arguments

- x:

  a \[rbi::libbi()\] object which has been run, or a list of data frames
  containing state trajectories (as returned by `bi_read`)

- origin:

  the time origin, i.e. the date or time corresponding to time 0

- unit:

  the unit of time that each time step corresponds to; this must be a
  unit understood by
  [`lubridate::period`](https://lubridate.tidyverse.org/reference/period.html),
  optionally with a number in advance, e.g. "day" or "2 weeks" or "3
  seconds"

- ...:

  any arguments for `bi_read` (e.g., `file`)

## Value

a list of data frames as returned by `bi_read`, but with real times
