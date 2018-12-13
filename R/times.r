##' Split a unit string
##'
##' Splits a unit string (e.g., "2 weeks") into the amount (2) and unit ("weeks")
##' @param unit_string the string to split
##' @return a list with two elements, "num" (the amount) and "unit", , for use with lubridate::period
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
  return(list(num=amount, unit=unit))
}

#' @rdname to_actual_times
#' @name to_actual_times
#' @title Convert LibBi times to actual times or dates
#' @description This function converts from LibBi times (i.e., 0, 1, 2, ...) to actual times or dates
#' @param x a \code{\link{libbi}} object which has been run, or a list of data
#'     frames containing state trajectories (as returned by \code{bi_read})
#' @param origin the time origin, i.e. the date or time corresponding to LibBi
#'   time 0
#' @param unit the unit of time that each time step corresponds to; this must bg
#'   a unit understood by \code{lubridate::period}, optionally with a number in
#'   advance, e.g. "day" or "2 weeks" or "3 seconds"
#' @param ... any arguments for \code{bi_read} (e.g., \code{file})
#' @return a list of data frames as returned bi \code{bi_read}, but with real times
#' @importFrom lubridate period
#' @importFrom rbi bi_read
#' @export
to_actual_times <- function(x, origin, unit, ...) {

  if (("libbi" %in% class(x)) || (is.character(x))) {
    vars <- do.call(bi_read, list(x=x, ...))
  } else if (is.list(x)) {
    vars <- x
  } else {
    stop("'x' must be a 'libbi' object or a file name or a list of data frames.")
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

  return(vars)
}

#' @rdname to_libbi_times
#' @name to_libbi_times
#' @title Convert actual times or dates to LibBi times
#' @description This function converts from real times/dates to LibBi times
#' @param x a data frame containing a "time" column, or a list containing such data frames
#' @param origin the time origin, i.e. the date or time corresponding to LibBi
#'   time 0
#' @param unit the unit of time that each time step corresponds to; this must bg
#'   a unit understood by \code{lubridate::period}, optionally with a number in
#'   advance, e.g. "day" or "2 weeks" or "3 seconds"
#' @return a list of data frames that can be passed to libbi
#' @importFrom lubridate as.interval period
#' @export
to_libbi_times <- function(x, origin, unit) {

  if (is.data.frame(x)) {
    vars <- list(df=x)
    list_given <- FALSE
  } else if (is.list(x)) {
    vars <- x
    list_given <- TRUE
  } else {
    stop("'x' must be a data frame or a list of data frames.")
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

  if (!list_given) { vars <- vars$df }

  return(vars)
}
