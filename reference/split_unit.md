# Split a unit string

Splits a unit string (e.g., "2 weeks") into the amount (2) and unit
("weeks")

## Usage

``` r
split_unit(unit_string)
```

## Arguments

- unit_string:

  the string to split

## Value

a list with two elements, "num" (the amount) and "unit", , for use with
lubridate::period

## Author

Sebastian Funk
