#' @rdname plot
#' @name plot
#' @title Plot results from libbi
#'
#' @description
#' Plots state trajectories (unless plot = FALSE) and invisibly returns a list of state trajectories and other plots.
#' @param x A \code{libbi} object containing Monte-Carlo samples
#' @param type character vector determining which plots to generate; options are: "state", "obs", "param", "noise", "logevals"; by default, all will be plotted; more specific selections of variables can be given as arguments with the name of the type containing character vectors of variables, e.g. \code{param="alpha"} to just plot parameter alpha (requiring "param" to be given as one of "type")
#' @param quantiles if plots are produced, which quantile to use for confidence intervals (NULL for no confidence intervals)
#' @param date.origin date of origin (if dates are to be calculated)
#' @param date.unit unit of date (if desired, otherwise the time dimension will be used); possible options: "day", "week", "month", "year", optionally with a quantity, e.g., "2 week"
#' @param data observations (a named list of data frames, a \code{libbi} object with observations, or a NetCDF file name)
#' @param extra.aes extra aesthetics (for ggplot)
#' @param all.times whether to plot all times (not only ones with observations)
#' @param hline horizontal marker lines, named vector in format (state = value)
#' @param burn How many iterations to burn
#' @param steps whether to plot lines as stepped lines
#' @param select list of selection criteria, as named list of selected elements. If the list contains "np", it is treated specially.
#' @param threshold thresholds for any of the trajectory variables; named list of named vector with a "lower" and/or "upper" threshold; any values outside this range for the variable with the name of the list element will not be plotted
#' @param data.colour colour for plotting the data
#' @param base.alpha base alpha value for credible intervals
#' @param np.alpha alpha of trajectories, if 'np' is part of \code{select} (default: 0.35)
#' @param trend how the trend should be characterised (e.g., mean, median, or NULL for no trend line)
#' @param limit.to.data whether to limit the time axis to times with observations (default: FALSE)
#' @param labels facet labels, in case they are to be rewritten, to be parsed using \code{label_parsed}; should be given as named character vector of (parameter = 'label') pairs
#' @param verbose if set to TRUE, additional output will be displayed
#' @param plot set to FALSE to suppress plot of trajectories
#' @param ... more specific selection of variables to plot (see the \code{type} option); any other options will be interpreted as options for geom_step / geom_line / geom_point / etc. when plotting states/noises/observations, e.g. lwd or others
#' @return a plot of trajectories
#' @import ggplot2 data.table
#' @importFrom lubridate %m+% years
#' @importFrom rbi bi_read bi_contents var_names
#' @importFrom stats quantile as.formula
#' @export
#' @examples
#' example_run_file <- system.file(package="rbi", "example_output.nc")
#' example_model_file <- system.file(package="rbi", "PZ.bi")
#' example_bi <- attach_data(libbi(example_model_file), "output", example_run_file)
#'
#' plot(example_bi) # just plot trajectories
#' \dontrun{p <- plot(example_bi, plot = FALSE) # get whole suite of plots
#'
#' p$trajectories
#' p$densities
#' p$traces
#' p$logevals}
#' @author Sebastian Funk
plot.libbi <- function(x, ..., 
                       type = c("state", "noise", "obs", "param", "logeval"),
                       quantiles = c(0.5, 0.95),
                       date.origin, date.unit,
                       data, extra.aes,
                       all.times = FALSE, hline,
                       burn, steps = FALSE, select, threshold,
                       data.colour = "red", base.alpha = 0.5,
                       trend = "median", np.alpha=0.35, limit.to.data = FALSE,
                       labels, verbose = FALSE, plot = TRUE)
{
    retval <- NULL
    ret_data <- list() ## list holding data to be returned

    dot_options <- list(...)

    use_dates <- FALSE
    summarise_columns <- c("np", "time", "time_next")

    all_types <- eval(formals(plot.libbi)[["type"]])
    missing_types <- setdiff(type, all_types)
    if (length(missing_types) > 0) {
        stop("Invalid 'type' argument(s): ", paste(missing_types) )
    }

    given_vars <- list()
    for (type.loop in intersect(names(dot_options), all_types)) {
        given_vars[[type.loop]] <- dot_options[[type.loop]]
        dot_options[[type.loop]] <- NULL
    }

    ## check plots: https://stat.ethz.ch/pipermail/r-help/2007-December/149117.html

    if (!missing(select) && "np" %in% names(select)) select_id <- TRUE
    else select_id <- FALSE

    if (missing(date.origin))
    {
        if (!missing(date.unit) && date.unit == "year")
        {
            use_dates <- TRUE
        }
    } else {
        if (!missing(date.unit))
        {
            use_dates <- TRUE
        } else {
            warning("date.origin given but no date.unit, will use time.dim instead")
        }
    }
    if (missing(labels)) labels <- c()

    data_missing <- missing(data)
    if (data_missing)
    {
        if (!is.null(x[["options"]]) &&
            !is.null(x[["options"]][["obs-file"]]))
        {
            ## if obs is missing but a libbi object passed, get obs file from
            ## the object
            vars_in_obs_file <- bi_contents(x[["options"]][["obs-file"]])
            data <- bi_read(x[["options"]][["obs-file"]],
                            vars=intersect(var_names(x[["model"]], type="obs"),
                                           vars_in_obs_file),
                            dims=x[["dims"]],
                            coord_dims=x[["coord_dims"]])
        }
    }

    if (steps)
    {
        ribbon_func <- function(mapping, ...)
        {
            if (missing(mapping)) mapping <- list()
            mapping <- utils::modifyList(mapping, aes_string(xmin = "time"))
            mapping <- utils::modifyList(mapping, aes_string(xmax = "time_next"))
            geom_rect(mapping, ...)
        }
        line_func <- geom_step
    } else
    {
        ribbon_func <- geom_ribbon
        line_func <- geom_line
    }

    if (!missing(extra.aes))
    {
        if ("color" %in% names(extra.aes) && !("fill" %in% names(extra.aes)))
        {
            extra.aes["fill"] <- extra.aes["color"]
        } else if ("fill" %in% names(extra.aes) & !("color" %in% names(extra.aes)))
        {
            extra.aes["color"] <- extra.aes["fill"]
        }
        summarise_columns <- c(summarise_columns, unique(extra.aes))
    }

    if (missing(burn)) burn <- 0

    existing_vars <- bi_contents(x)

    vars <- list()
    init_vars <- c()
    for (type.loop in type)
    {
        if (!(type.loop %in% names(given_vars)))
        {
            if (type.loop == "logeval")
            {
                given_vars[[type.loop]] <- intersect(existing_vars, c("loglikelihood", "logprior", "logweight", "logevidence"))
            } else
            {
                type_vars <- var_names(x[["model"]], type=type.loop)
                if (type.loop == "param")
                {
                    init_vars <-
                        sub("[[:space:]]*~.*$", "",
                            grep("~", get_block(x[["model"]], "proposal_initial"),
                                 value=TRUE))
                    init_vars <- setdiff(init_vars, paste0(type_vars, "_0"))
                    type_vars <- c(type_vars, init_vars)
                }
                given_vars[[type.loop]] <-
                    intersect(existing_vars, type_vars)
            }
        }
        vars[[type.loop]] <- intersect(given_vars[[type.loop]], existing_vars)
        missing_vars <- setdiff(given_vars[[type.loop]], vars[[type.loop]])
        if (length(missing_vars) > 0)
        {
            warning("Variable(s) ", missing_vars, " not found in given 'x'.")
        }
    }

    ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@7@"]]));##:ess-bp-end:##
    missing_types <- setdiff(names(given_vars), type)
    if (length(missing_types) > 0) {
        warning("Variables given for type(s) ", paste(missing_types), ", but not included in 'type' variable. Will not plot these")
    }

    trajectory_vars <- c("state", "obs", "noise")
    vars[["trajectories"]] <- unname(unlist(vars[trajectory_vars]))
    vars[trajectory_vars] <- NULL

    clean_data <- function(y, name, file_name, ...)
    {
        if (missing(file_name))
        {
            y <- rbi::bi_read(y, ...)
        } else
        {
            opt_name <- paste(file_name, "file", sep="-")
            if (!is.null(y[["options"]]) &&
                !is.null(y[["options"]][[opt_name]]))
            {
                y <- rbi::bi_read(y, ...)
            } else
            {
                stop("No ", opt_name, " found.")
            }
        }
        y <- lapply(y, function(z) {
          if (is.data.frame(z)) { data.table::data.table(z) } else {z}
        })
    }

    clean_dates <- function(values, time.dim, use_dates, date.unit, date.origin)
    {
        if (time.dim %in% colnames(values))
        {
            if (use_dates)
            {
                known_units <- c("day", "week", "month", "year")
                date_split <- unlist(strsplit(date.unit, " "))
                if (length(date_split) == 1)
                {
                    amount <- 1
                    unit <- date_split
                } else if (length(date_split) == 2)
                {
                    amount <- as.numeric(date_split[1])
                    unit <- date_split[2]
                } else
                {
                    stop("'date.unit' must consist of one or two words")
                }
                if (grepl("days?$", unit))
                {
                    values[, paste("time") := date.origin + get(time.dim) * amount]
                    values[, paste("time_next") := get("time") + 1 * amount]
                } else if (grepl("weeks?$", unit))
                {
                    values[, paste("time") := date.origin + get(time.dim) * 7 * amount]
                    values[, paste("time_next") := get("time") + 7 * amount]
                } else if (grepl("months?$", unit))
                {
                    values[, paste("time") := date.origin %m+% months(as.integer(get(time.dim) * amount))]
                    values[, paste("time_next") := get("time") %m+% months(amount)]
                } else if (grepl("years?$", unit))
                {
                    if (missing(date.origin)) {
                        values[, paste("time") := as.Date(paste(get(time.dim) * amount, 1, 1, sep = "-"))]
                        values[, paste("time_next") := as.Date(paste(get(time.dim) + amount, 1, 1, sep = "-"))]
                    } else {
                        values[, paste("time") := date.origin + years(as.integer(get(time.dim) * amount))]
                        values[, paste("time_next") := get("time") %m+% months(12 * amount)]
                    }
                } else {
                    stop("Unknown date unit: ", unit)
                }
            } else {
                values[, paste("time") := get(time.dim)]
                values[, paste("time_next") := get("time") + 1]
            }
        }
        return(values)
    }

    time_dim <- ifelse(length(x$time_dim) == 1, x$time_dim, "time")

    if (!missing(select) && time_dim %in% names(select)) {
      temp_time_df <- data.table(time=select[[time_dim]])
      temp_time_df <- clean_dates(temp_time_df, time_dim, use_dates, date.unit, date.origin)
      select[[time_dim]] <- temp_time_df[[time_dim]]
    }

    ## plot trajectories
    remove_times <- list()
    if (length(intersect(type, c("state", "obs", "noise"))) > 0)
    {
        vdt <- NULL
        if (length(vars[["trajectories"]]) > 0) {
            if (verbose) message(date(), " Getting trajectory samples")
            samples <- clean_data(x, "x", verbose=verbose, vars=vars[["trajectories"]])
        }
        if (!data_missing) {
          if (verbose) message(date(), " Getting observations")
            data <- clean_data(data, "data", "obs", verbose=verbose)
        }

        if (verbose) message(date(), " Generating trajectory plots")
        for (var in vars[["trajectories"]])
        {
            if (var %in% names(samples) && nrow(samples[[var]]) > 0)
            {
                values <- samples[[var]]
                if ("np" %in% colnames(samples[[var]]) && burn > 0)
                {
                    values <- values[get("np") >= burn]
                    if (nrow(values) == 0) {
                        stop("Nothing left after burn-in")
                    }
                }

                values <- clean_dates(values, time_dim, use_dates, date.unit, date.origin)

                if (!missing(select))
                {
                    for (var_name in names(select))
                    {
                        if (var_name %in% colnames(values))
                        {
                            if (!(var_name == "np") || time_dim %in% colnames(values)) {
                                values <- values[get(var_name) %in% select[[var_name]]]
                                if (class(values[, get(var_name)]) == "factor")
                                {
                                    values[, paste(var_name) := factor(get(var_name))]
                                }
                            }
                        }
                    }
                }

                if (!missing(threshold) && (var %in% names(threshold)))
                {
                    if ("lower" %in% names(threshold[[var]]))
                    {
                        lower <- threshold[[var]][["lower"]]
                        remove_times[[var]] <- values[get("value") < lower]
                    }
                    if ("upper" %in% names(threshold[[var]]))
                    {
                        upper <- threshold[[var]][["upper"]]
                        remove_times[[var]] <- values[get("value") > upper]
                    }
                }

                sum.by <- intersect(summarise_columns, colnames(values))
                values <- values[, list(value = sum(get("value"))), by = sum.by]

                state.wo <- setdiff(setdiff(summarise_columns, "np"),
                                    colnames(values))
                for (wo in state.wo)
                {
                    values[, paste(wo) := "n/a"]
                }

                new_vars <- data.table::data.table(var = rep(var, nrow(values)), values)
                if (is.null(vdt))
                {
                    vdt <- new_vars
                } else
                {
                    vdt <- rbind(vdt, new_vars)
                }
            } else if (var != ".var")
            {
                warning(paste("Variable", var, "does not exist"))
            }
        }

        ## factorise columns
        if (!is.null(vdt)) vdt <- factorise_columns(vdt, labels)

        if (!is.null(vdt) && nrow(vdt) > 0 && length(remove_times) > 0)
        {
            remove_times <- rbindlist(remove_times)
            remove_times <- unique(remove_times)
            remove_times[, paste("value") := NULL]
            remove_times[, paste("remove") := TRUE]
            vdt <- merge(vdt, remove_times, all.x=TRUE)
            vdt <- vdt[is.na(get("remove"))]
            vdt[, paste("remove") := NULL]
        }

        if (!missing(data))
        {
            dataset <- lapply(names(data), function(y) {data.table::data.table(data[[y]])[, var := y]})
            dataset <- rbindlist(dataset, fill=TRUE)
            dataset <- factorise_columns(dataset, labels)
            dataset <- clean_dates(dataset, time_dim, use_dates, date.unit, date.origin)

            if (!all.times && !is.null(vdt) && nrow(vdt) > 0)
            {
                if (limit.to.data)
                {
                    ## for all states, only retain times with observations
                    vdt <- vdt[get("time") %in% dataset[, get("time")]]
                } else
                {
                    for (data_var in unique(dataset[, var]))
                    {
                        ## for states in observations, only retain times with observations
                        vdt <- vdt[(var != data_var) |
                                   (get("time") %in% dataset[var == data_var, get("time")])]
                    }
                }
            }

            dataset <- dataset[var %in% vars[["trajectories"]]]
            if (nrow(dataset) > 0)
            {
                for (col in colnames(dataset))
                {
                    dataset[is.na(get(col)), paste(col) := "n/a"]
                }
                if (!missing(select))
                {
                    for (var_name in names(select))
                    {
                        if (var_name %in% colnames(dataset))
                        {
                            dataset <- dataset[get(var_name) %in% select[[var_name]]]
                            if (class(values[, get(var_name)]) == "factor")
                            {
                                dataset[, paste(var_name) := factor(get(var_name))]
                            }
                        }
                    }
                }

                for (i in seq_along(quantiles))
                {
                    dataset[, paste("min", i, sep = ".") := 0]
                    dataset[, paste("max", i, sep = ".") := 0]
                }

                for (missing_column in setdiff(colnames(vdt), colnames(dataset)))
                {
                    dataset[, paste(missing_column) := "n/a"]
                }
            }
            ret_data[["observations"]] <- dataset
        }

        aggregate_values <- NULL
        if (!is.null(vdt) && nrow(vdt) > 0) {
            var.by <- c("var", intersect(setdiff(summarise_columns, "np"), colnames(vdt)))

            if (!is.null(trend))
            {
                aggregate_values <- vdt[, list(value = do.call(trend, list(get("value"), na.rm = TRUE))), by = var.by]
            }

            if (!is.null(quantiles))
            {
                for (i in seq_along(quantiles))
                {
                    quantile_values <-
                        vdt[, list(max = stats::quantile(get("value"), 0.5 + quantiles[i] / 2, na.rm = TRUE),
                                   min = stats::quantile(get("value"), 0.5 - quantiles[i] / 2, na.rm = TRUE)),
                            by = var.by]
                    data.table::setnames(quantile_values, c("min", "max"), paste(c("min",  "max"),  i,  sep = "."))
                    if (is.null(aggregate_values))
                    {
                        aggregate_values <- quantile_values
                    } else
                    {
                        aggregate_values <- merge(aggregate_values, quantile_values, by = var.by)
                    }
                }
                if (steps)
                {
                    max.time <- aggregate_values[, max(get("time"))]
                    for (i in seq_along(quantiles))
                    {
                        aggregate_values[get("time") == max.time,
                                         paste("min", i, sep = ".") := NA]
                        aggregate_values[get("time") == max.time,
                                         paste("max", i, sep = ".") := NA]
                    }
                }
            }
            if (!is.null(aggregate_values))
            {
                vars_n <- aggregate_values[, list(single = (.N == 1)), by = var]
                aggregate_values <- merge(aggregate_values, vars_n, by = "var", all.x = TRUE)
                ret_data[["trajectories"]] <- aggregate_values[, !"single", with = FALSE]
            }
        }

        if (!is.null(vdt) && nrow(vdt) > 0)
        {
          vars_n <- vdt[, list(single = (.N == 1)), by = var]
          vdt <- merge(vdt, vars_n, by = "var", all.x = TRUE)
        }

        aesthetic <- list(x = "time", y = "value")
        if (!missing(extra.aes))
        {
            aesthetic <- c(aesthetic, extra.aes)
        }

        if (!missing(extra.aes) && "color" %in% names(extra.aes) && select_id && !is.null(vdt) && nrow(vdt) > 0)
        {
            vdt <- vdt[, paste("color_np") :=
                             paste(get(extra.aes["color"]),
                                   get("np"), sep = "_")]
        }

        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@6@"]]));##:ess-bp-end:##
        p <- ggplot(mapping = do.call(aes_string, aesthetic))

        if (!is.null(quantiles) && !is.null(aggregate_values) && nrow(aggregate_values) > 0)
        {
            alpha <- base.alpha
            for (i in seq_along(quantiles))
            {
                str <- as.list(paste(c("max", "min"), i, sep = "."))
                names(str) <- c("ymax", "ymin")
                p <- p + ribbon_func(data = aggregate_values[single == FALSE], do.call(aes_string, str), alpha = alpha)
                alpha <- alpha / 2
                if (!is.null(aggregate_values) && nrow(aggregate_values[single == TRUE]) > 0 && i == 1) ## only plot first as errorbar
                {
                    p <- p + do.call(geom_errorbar, c(list(data = aggregate_values[single == TRUE], do.call(aes_string, str)), dot_options))
                }
            }
        }
        if ("color" %in% names(aesthetic))
        {
            if (!is.null(trend) && !is.null(aggregate_values) && nrow(aggregate_values) > 0)
            {
                p <- p + do.call(line_func, c(list(data = aggregate_values[single == FALSE]), dot_options))
                if (nrow(aggregate_values[single == TRUE]) > 0)
                {
                    p <- p + do.call(geom_point, c(list(data = aggregate_values[single == TRUE], shape = 4), dot_options))
                }
            }
            if (select_id && !is.null(vdt) && nrow(vdt) > 0)
            {
                p <- p + do.call(line_func, c(list(data = vdt[single == FALSE], mapping = aes_(group =~ color_np), alpha = np.alpha), dot_options))
                if (nrow(vdt[single == TRUE]) > 0)
                {
                    p <- p + do.call(geom_point, c(list(data = vdt[single == TRUE], aes_(group =~ factor(np)), shape = 4, alpha = np.alpha), dot_options))
                }
            }
        } else
        {
            if (!is.null(trend) && !is.null(aggregate_values) && nrow(aggregate_values) > 0)
            {
                p <- p + do.call(line_func, c(list(data = aggregate_values[single == FALSE]), dot_options))
                p <- p + do.call(geom_point, c(list(data = aggregate_values[single == TRUE], shape = 4), dot_options))
            }
            if (select_id && !is.null(vdt) && nrow(vdt) > 0)
            {
                p <- p + do.call(line_func, c(list(data = vdt[single == FALSE], aes_(group =~ factor(np)), alpha = np.alpha), dot_options))
                p <- p + do.call(geom_point, c(list(data = vdt[single == TRUE], aes_(group =~ factor(np)), alpha = np.alpha, shape = 4), dot_options))
            }
        }
        p <- p + scale_y_continuous() + ylab("")
        ## p <- p + expand_limits(y = 0)
        if (!missing(data) && nrow(dataset) > 0)
        {
            if (!missing(extra.aes) && "color" %in% names(extra.aes))
            {
                p <- p + geom_point(data = dataset)
            } else
            {
                p <- p + geom_point(data = dataset, color = data.colour, size = 2)
            }
        }
        p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                       legend.position = "top")
        if (use_dates)
        {
            p <- p + scale_x_date("")
        }
        if (!missing(hline))
        {
          named <- which(names(hline) != "")
          if (is.null(names(hline)))
          {
            unnamed <- seq_along(hline)
          } else
          {
            unnamed <- which(names(hline) == "")
          }
          for (hline_var_id in named)
          {
            if (names(hline)[hline_var_id] %in% vars$trajectories) {
              hline_data <- data.frame(var = names(hline)[hline_var_id],
                                       yintercept = hline[hline_var_id])
              p <- p + geom_hline(data = hline_data,
                                  aes_(yintercept =~ yintercept), color = "black")
            }
          }
          for (hline_var_id in unnamed)
          {
            hline_data <- data.frame(yintercept = hline[hline_var_id])
            p <- p + geom_hline(data = hline_data,
                                aes_(yintercept =~ yintercept), color = "black")
          }
        }

        if (length(vars[["trajectories"]]) > 1 &&
            ((!is.null(aggregate_values) && nrow(aggregate_values) > 0) ||
             (!is.null(vdt) && nrow(vdt) > 0)))
        {
          p <- p +
            facet_wrap(~ var, scales = "free_y",
                       ncol = round(sqrt(length(vars[["trajectories"]]))), 
                       labeller = label_parsed)
        }

        retval <- p
    }

    ## return all plots invisibly
    invisible(retval)
}

