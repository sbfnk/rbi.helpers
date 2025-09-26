##' Split a unit string
##'
##' Splits a unit string (e.g., "2 weeks") into the amount (2) and unit
##' ("weeks")
##' @param unit_string the string to split
##' @return a list with two elements, "num" (the amount) and "unit", , for use
##'   with lubridate::period
##' @author Sebastian Funk
##' @keywords internal
split_unit <- function(unit_string) {
  unit_string <- sub("^[[:space:]]*", "", unit_string)
  unit_string <- sub("[[:space:]]*$", "", unit_string)
  ## check if given as number + unit
  if (grepl("^[0-9]+[[:space:]]*", unit_string)) {
    split_unit <- strsplit(unit_string, "[[:space:]]+")[[1]]
    amount <- as.numeric(split_unit[1])
    unit <- split_unit[2]
  } else {
    amount <- 1
    unit <- unit_string
  }
  list(num = amount, unit = unit)
}

#' @rdname numeric_to_time
#' @name numeric_to_time
#' @title Convert numeric times to actual times or dates
#' @description This function converts from numeric times (i.e., 0, 1, 2, ...)
#'   to actual times or dates
#' @param x a [rbi::libbi()] object which has been run, or a list of data
#'     frames containing state trajectories (as returned by \code{bi_read})
#' @param origin the time origin, i.e. the date or time corresponding to time 0
#' @param unit the unit of time that each time step corresponds to; this must be
#'   a unit understood by \code{lubridate::period}, optionally with a number in
#'   advance, e.g. "day" or "2 weeks" or "3 seconds"
#' @param ... any arguments for \code{bi_read} (e.g., \code{file})
#' @return a list of data frames as returned by \code{bi_read}, but with real
#'   times
#' @importFrom lubridate period
#' @importFrom rbi bi_read
#' @export
numeric_to_time <- function(x, origin, unit, ...) {
  if (("libbi" %in% class(x)) || (is.character(x))) {
    vars <- do.call(bi_read, list(x = x, ...))
  } else if (is.list(x)) {
    vars <- x
  } else {
    stop(
      "'x' must be a 'libbi' object or a file name or a list of data frames."
    )
  }

  ## convert unit string to time step for lubridate::period
  time_step <- split_unit(unit)

  for (var in names(vars)) {
    ## check if data frame has a time variable
    if (is.data.frame(vars[[var]]) && "time" %in% colnames(vars[[var]])) {
      vars[[var]][["time"]] <-
        origin + vars[[var]][["time"]] * do.call(lubridate::period, time_step)
    }
  }

  vars
}

#' @rdname time_to_numeric
#' @name time_to_numeric
#' @title Convert actual times or dates to numeric times
#' @description This function converts from real times/dates to numeric times
#'   (0, 1, 2, ...)
#' @param x a data frame containing a "time" column, or a list containing such
#'   data frames
#' @param origin the time origin, i.e. the date or time corresponding to time 0
#' @param unit the unit of time that each time step corresponds to; this must be
#'   a unit understood by \code{lubridate::period}, optionally with a number in
#'   advance, e.g. "day" or "2 weeks" or "3 seconds"
#' @return a list of data frames that can be passed to \code{libbi}
#' @importFrom lubridate as.interval period
#' @export
time_to_numeric <- function(x, origin, unit) {
  if (is.data.frame(x)) {
    vars <- list(df = x)
    list_given <- FALSE
  } else if (is.list(x)) {
    vars <- x
    list_given <- TRUE
  } else {
    stop("'x' must be a data frame or a list.")
  }

  ## convert unit string to time step for lubridate::period
  time_step <- split_unit(unit)

  for (var in names(vars)) {
    ## check if data frame has a time variable
    if (is.data.frame(vars[[var]]) && "time" %in% colnames(vars[[var]])) {
      vars[[var]][["time"]] <-
        as.interval(vars[[var]][["time"]] - origin, origin) /
        do.call(period, time_step)
    }
  }

  if (!list_given) {
    vars <- vars$df
  }

  vars
}
