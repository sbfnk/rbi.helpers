# Convert actual times or dates to numeric times

This function converts from real times/dates to numeric times (0, 1, 2,
...)

## Usage

``` r
time_to_numeric(x, origin, unit)
```

## Arguments

- x:

  a data frame containing a "time" column, or a list containing such
  data frames

- origin:

  the time origin, i.e. the date or time corresponding to time 0

- unit:

  the unit of time that each time step corresponds to; this must be a
  unit understood by
  [`lubridate::period`](https://lubridate.tidyverse.org/reference/period.html),
  optionally with a number in advance, e.g. "day" or "2 weeks" or "3
  seconds"

## Value

a list of data frames that can be passed to `libbi`
