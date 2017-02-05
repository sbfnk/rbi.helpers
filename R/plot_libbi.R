#' Plot results from libbi
#'
#' Plots state trajectories (unless plot = FALSE) and invisibly returns a list of state trajectories and other plots.
#' @param x Monte-Carlo samples, either a \code{libbi} object or a list of data frames, as returned by \code{bi_read}, or the name of an NetCDF file used as 'output' in a libbi run
#' @param model model file or a \code{bi_model} object (if \code{data} is not a \code{libbi} object)
#' @param prior optional; Prior samples, either a \code{libbi} object or a list of data frames, as returned by \code{bi_read}
#' @param type character vector determining which plots to generate; options are: "state", "obs", "param", "noise", "logevals"; by default, all will be plotted; more specific selections of variables can be given as arguments with the name of the type containing character vectors of variables, e.g. \code{param="alpha"} to just plot parameter alpha (requiring "param" to be given as one of "type")
#' @param logevals logged density evaluations (e.g., loglikelihood, logprior, logweights, logevidence) to plot; if not given, all logged density evaluations will be plotted
#' @param quantiles if plots are produced, which quantile to use for confidence intervals (NULL for no confidence intervals)
#' @param date.origin date of origin (if dates are to be calculated)
#' @param date.unit unit of date (if desired, otherwise the time dimension will be used)
#' @param time.dim time dimension ("time" by default)
#' @param data observations (a named list of data frames, a \code{libbi} object with observations, or a NetCDF file name)
#' @param extra.aes extra aesthetics (for ggplot)
#' @param all.times whether to plot all times (not only ones with observations)
#' @param hline horizontal marker lines, named vector in format (state = value)
#' @param burn How many iterations to burn
#' @param steps whether to plot lines as stepped lines
#' @param select list of selection criteria, as named list of selected elements. If the list contains "np", it is treated specially.
#' @param data.colour colour for plotting the data
#' @param base.alpha base alpha value for credible intervals
#' @param np.alpha alpha of trajectories, if 'np' is part of \code{select} (default: 0.35)
#' @param trend how the trend should be characterised (e.g., mean, median, or NULL for no trend line)
#' @param densities density geometry (e.g., "histogram" (default) or "density")
#' @param density_args list of arguments to pass to density geometry
#' @param limit.to.data whether to limit the time axis to times with observations
#' @param labels facet labels, in case they are to be rewritten, to be parsed using \code{label_parsed}; should be given as named character vector of (parameter = 'label') pairs
#' @param brewer.palette optional; brewer color palette
#' @param verbose if set to TRUE, additional output will be displayed
#' @param plot set to FALSE to suppress plot of trajectories
#' @param ... more specific selection of variables to plot (see the \code{type} option); any other options will be interpreted as options for geom_step / geom_line / geom_point / etc. when plotting states/noises/observations, e.g. lwd or others 
#' @return a list of plots: states, densities, traces, correlations, noises, logdensities, as well as underlying raw and aggregate data
#' @import ggplot2 scales data.table
#' @importFrom lubridate %m+% years
#' @importFrom rbi bi_read bi_model
#' @importFrom stats quantile as.formula
#' @importFrom GGally ggcorr
#' @export
#' @examples
#' example_run_file <- system.file(package="rbi.helpers", "example_run.nc")
#' example_model_file <- system.file(package="rbi", "PZ.bi")
#' example_bi <- recreate_libbi(example_model_file, example_run_file)
#'
#' plot(example_bi) # just plot trajectories
#' p <- plot(example_bi, plot = FALSE) # get whole suite of plots
#'
#' p$trajectories
#' p$correlations
#' p$pairs
#' p$densities
#' p$traces
#' p$logevals
#' @author Sebastian Funk
plot_libbi <- function(x, model, prior,
                       type = c("state", "noise", "obs", "param", "logeval"),
                       quantiles = c(0.5, 0.95),
                       date.origin, date.unit, time.dim = "time",
                       data, extra.aes,
                       all.times = FALSE, hline,
                       burn, steps = FALSE, select,
                       data.colour = "red", base.alpha = 0.5,
                       np.alpha=0.35, trend = "median",
                       densities = "histogram",
                       density_args = list(), limit.to.data = FALSE,
                       labels, brewer.palette, verbose = FALSE,
                       plot = TRUE, ...)
{
    plots <- list() ## list holding the plots to be returned
    ret_data <- list() ## list holding data to be returned

    dot_options <- list(...)

    use_dates <- FALSE
    summarise_columns <- c("np", "time", "time_next")

    all_types <- eval(formals(plot_libbi)[["type"]])
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

    if (missing(data))
    {
        if ("libbi" %in% class(x) &&
            !is.null(x[["options"]]) &&
            !is.null(x[["options"]][["obs-file"]]))
        {
            ## if obs is missing but a libbi object passed, get obs file from
            ## the object
            vars_in_obs_file <- bi_contents(x[["options"]][["obs-file"]])
            data <- bi_read(x[["options"]][["obs-file"]],
                            vars=intersect(var_names(x[["model"]], "obs"),
                                           vars_in_obs_file))
        }
    }

    if (missing(model))
    {
        if ("libbi" %in% class(x))
        {
            model <- x$model
        }
    } else
    {
        if ("libbi" %in% class(x))
        {
            stop("'model' should not be given if 'x' is a 'libbi' object'.")
        }
        if (is.character(model)) {
            model <- rbi::bi_model(model)
        } else if (!("bi_model" %in% class(model))) {
            stop("'model' must be either a 'bi_model' object or a path to a valid model file in LibBi's syntax")
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

   if ("libbi" %in% class(x) || "character" %in% class(x)) {
      existing_vars <- bi_contents(x)
    } else {
      existing_vars <- names(x)
    }

    vars <- list()
    for (type.loop in type)
    {
        if (!(type.loop %in% names(given_vars)))
        {
            if (missing(model))
            {
                given_vars[[type.loop]] <- c()
            } else
            {
                if (type.loop == "logeval")
                {
                    given_vars[[type.loop]] <- intersect(existing_vars, c("loglikelihood", "logprior", "logweight", "logevidence"))
                } else
                {
                    given_vars[[type.loop]] <- intersect(existing_vars, var_names(model, type.loop))
                }
            }
        }
        vars[[type.loop]] <- intersect(given_vars[[type.loop]], existing_vars)
        missing_vars <- setdiff(given_vars[[type.loop]], vars[[type.loop]])
        if (length(missing_vars) > 0)
        {
            warning("Variable(s) ", missing_vars, " not found in given 'x'.")
        }
    }

    missing_types <- setdiff(names(given_vars), type)
    if (length(missing_types) > 0) {
        warning("Variables given for type(s) ", paste(missing_types), ", but not included in 'type' variable. Will not plot these")
    }

    trajectory_vars <- c("state", "obs", "noise")
    vars[["trajectories"]] <- unname(unlist(vars[trajectory_vars]))
    vars[trajectory_vars] <- NULL

    clean_data <- function(y, name, file_name, ...)
    {
        if ("libbi" %in% class(y))
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
        } else if (is.data.frame(y))
        {
            y <- list(.var = y)
        } else if (is.character(y))
        {
            y <- rbi::bi_read(y, ...)
        } else if (!is.list(y))
        {
            stop("'", name, "' must be a 'libbi' object or a list of data frames or a data frame.")
        }
        y <- lapply(y, function(z) { if (is.data.frame(z)) { data.table::data.table(z) } else {z} })
    }

    clean_dates <- function(values, time.dim, use_dates, date.unit, date.origin)
    {
        if (time.dim %in% colnames(values))
        {
            if (use_dates)
            {
                if (date.unit == "day")
                {
                    values[, time := date.origin + get(time.dim)]
                    values[, time_next := time + 1]
                } else if (date.unit == "week")
                {
                    values[, time := date.origin + get(time.dim) * 7]
                    values[, time_next := time + 7]
                } else if (date.unit == "month")
                {
                    values[, time := date.origin %m+% months(as.integer(get(time.dim)))]
                    values[, time_next := time %m+% months(1)]
                } else if (date.unit == "year")
                {
                    if (missing(date.origin)) {
                        values[, time := as.Date(paste(get(time.dim), 1, 1, sep = "-"))]
                        values[, time_next := as.Date(paste(get(time.dim) + 1, 1, 1, sep = "-"))]
                    } else {
                        values[, time := date.origin + years(as.integer(get(time.dim)))]
                        values[, time_next := time %m+% months(12)]
                    }
                }
            } else {
                values[, time := get(time.dim)]
                values[, time_next := time + 1]
            }
        }
        return(values)
    }

    ## plot trajectories
    if (length(intersect(type, c("state", "obs", "noise"))) > 0)
    {
        vdt <- NULL
        if (length(vars[["trajectories"]]) > 0) {
            if (verbose) message(date(), " Getting trajectory samples")
            samples <- clean_data(x, "x", verbose=verbose, vars=vars[["trajectories"]])
        }
        if (!missing(prior)) {
            if (verbose) message(date(), " Getting prior samples")
            prior <- clean_data(prior, "prior", verbose=verbose)
        }
        if (!missing(data)) {
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
                    values <- values[np >= burn]
                    if (nrow(values) == 0) {
                        stop("Nothing left after burn-in")
                    }
                }

                if (!missing(select))
                {
                    for (var_name in names(select))
                    {
                        if (var_name %in% colnames(values))
                        {
                            if (!(var_name == "np") || time.dim %in% colnames(values)) {
                                values <- values[get(var_name) %in% select[[var_name]]]
                                if (class(values[, get(var_name)]) == "factor")
                                {
                                    values[, paste(var_name) := factor(get(var_name))]
                                }
                            }
                        }
                    }
                }

                values <- clean_dates(values, time.dim, use_dates, date.unit, date.origin)
                sum.by <- intersect(summarise_columns, colnames(values))
                values <- values[, list(value = sum(value)), by = sum.by]

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

        if (!missing(data))
        {
            data <- data[intersect(vars[["trajectories"]], names(data))]
            dataset <- lapply(names(data), function(y) {data.table::data.table(data[[y]])[, var := y]})
            dataset <- rbindlist(dataset, fill=TRUE)
            dataset <- factorise_columns(dataset, labels)
            dataset <- clean_dates(dataset, time.dim, use_dates, date.unit, date.origin)
            for (col in colnames(dataset))
            {
                dataset[is.na(get(col)), paste(col) := "n/a"]
            }

            if (nrow(dataset) > 0)
            {
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

                if (!all.times && !is.null(vdt) && nrow(vdt) > 0)
                {
                    if (limit.to.data)
                    {
                        ## for all states, only retain times with observations
                        vdt <- vdt[time %in% dataset[, time]]
                    } else
                    {
                        for (data_var in unique(dataset[, var]))
                        {
                            ## for states in observations, only retain times with observations
                            vdt <- vdt[(var != data_var) |
                                       (time %in% dataset[var == data_var, time])]
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
        }

        aggregate_values <- NULL
        if (!is.null(vdt) && nrow(vdt) > 0) {
            var.by <- c("var", intersect(setdiff(summarise_columns, "np"), colnames(vdt)))

            if (!is.null(trend))
            {
                aggregate_values <- vdt[, list(value = do.call(trend, list(value, na.rm = TRUE))), by = var.by]
            }

            if (!is.null(quantiles))
            {
                for (i in seq_along(quantiles))
                {
                    quantile_values <-
                        vdt[, list(max = stats::quantile(value, 0.5 + quantiles[i] / 2, na.rm = TRUE),
                                   min = stats::quantile(value, 0.5 - quantiles[i] / 2, na.rm = TRUE)),
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
                    max.time <- aggregate_values[, max(time)]
                    for (i in seq_along(quantiles))
                    {
                        aggregate_values[time == max.time, paste("min", i, sep = ".") := NA]
                        aggregate_values[time == max.time, paste("max", i, sep = ".") := NA]
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
          vdt <- vdt[, color_np := paste(get(extra.aes["color"]), get("np"), sep = "_")]
        }

        p <- ggplot(mapping = do.call(aes_string, aesthetic))

        if (!is.null(quantiles) && !is.null(aggregate_values) && nrow(aggregate_values) > 0)
        {
            alpha <- base.alpha
            for (i in seq_along(quantiles))
            {
                str <- as.list(paste(c("max", "min"), i, sep = "."))
                names(str) <- c("ymax", "ymin")
                p <- p + ribbon_func(data = aggregate_values, do.call(aes_string, str), alpha = alpha)
                alpha <- alpha / 2
                if (!is.null(aggregate_values) && nrow(aggregate_values[single == TRUE]) > 0)
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
                p <- p + do.call(line_func, c(list(data = vdt[single == FALSE], mapping = aes(group = color_np), alpha = np.alpha), dot_options))
                if (nrow(vdt[single == TRUE]) > 0)
                {
                    p <- p + do.call(geom_point, c(list(data = vdt[single == TRUE], aes(group = factor(np)), shape = 4, alpha = np.alpha), dot_options))
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
                p <- p + do.call(line_func, c(list(data = vdt[single == FALSE], aes(group = factor(np)), alpha = np.alpha), dot_options))
                p <- p + do.call(geom_point, c(list(data = vdt[single == TRUE], aes(group = factor(np)), alpha = np.alpha, shape = 4), dot_options))
            }
        }
        p <- p + scale_y_continuous(labels = comma) + ylab("")
        if (!missing(brewer.palette))
        {
            p <- p + scale_color_brewer(palette = brewer.palette)
            p <- p + scale_fill_brewer(palette = brewer.palette)
        }
        p <- p + expand_limits(y = 0)
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
            hline_data <- data.frame(var = names(hline)[hline_var_id],
                                     yintercept = hline[hline_var_id])
            p <- p + geom_hline(data = hline_data,
                                aes(yintercept = yintercept), color = "black")
          }
          for (hline_var_id in unnamed)
          {
            hline_data <- data.frame(yintercept = hline[hline_var_id])
            p <- p + geom_hline(data = hline_data,
                                aes(yintercept = yintercept), color = "black")
          }
        }

        if (length(vars[["trajectories"]]) > 1)
        {
          p <- p + facet_wrap(~ var, scales = "free_y",
                              ncol = round(sqrt(length(unlist(vars)))),
                              labeller = label_parsed)
        }

        plots[["trajectories"]] <- p
    }

    if ("param" %in% type)
    {
        pdt <- NULL
        if (length(vars[["param"]]) > 0) {
            if (verbose) message(date(), " Getting parameters")
            samples <- clean_data(x, "x", verbose=verbose, vars=vars[["param"]])
        }

        if (verbose) message(date(), " Parameter plots")
        for (param in vars[["param"]])
        {
            param_values <- list()
            param_values[["posterior"]] <- samples[[param]]
            if ("np" %in% colnames(samples[[param]]))
            {
                param_values[["posterior"]] <-
                    param_values[["posterior"]][np >= burn]
            }

            if (!missing(prior) && param %in% names(prior))
            {
                param_values[["prior"]] <- prior[[param]]
            }

            for (dist in names(param_values))
            {
                values <- param_values[[dist]]
                if (!("data.frame" %in% class(values)))
                {
                    values <- data.table::data.table(np = 0, value = values)
                }

                by.mean <- "np"
                if (!missing(extra.aes))
                {
                    by.mean <- c(by.mean, unique(unname(extra.aes)))
                }
                param.by <- intersect(by.mean, colnames(values))
                param.wo <- setdiff(by.mean, colnames(values))

                values <- values[, list(value = value), by = param.by]

                if (!("np" %in% colnames(values)))
                {
                    values[, np := 1]
                }

                for (wo in setdiff(param.wo, "np"))
                {
                    values[, paste(wo) := "n/a"]
                }

                new_params <- data.table::data.table(distribution = dist, parameter = rep(param, nrow(values)), values)
                if (is.null(pdt))
                {
                    pdt <- new_params
                } else
                {
                    pdt <- rbind(pdt, new_params)
                }
            }
        }
        if (!is.null(pdt) && nrow(pdt) > 0)
        {
            ## factorise columns
            pdt <- factorise_columns(pdt, labels)

            by.varying <- c("parameter", "distribution")
            if (!missing(extra.aes))
            {
                by.varying <- c(by.varying, unlist(unique(unname(extra.aes))))
            }
            pdt[, varying := (length(unique(value)) > 1), by = by.varying]

            ret_data <- c(ret_data, list(params = pdt))

            aesthetic <- list(x = "value", y = "..density..")
            if (!missing(extra.aes))
            {
                aesthetic <- c(aesthetic, extra.aes)
            }

            black_prior <- FALSE
            if (!missing(prior))
            {
                if (!missing(extra.aes) && length(intersect(c("color", "fill"), names(extra.aes))) > 0)
                {
                    ## if posterior is colourful, make prior black
                    black_prior <- TRUE
                } else
                {
                    aesthetic <- c(aesthetic, list(color = "distribution", fill = "distribution"))
                }
            }

            if (!is.null(pdt) && nrow(pdt[varying == TRUE & distribution == "posterior"]) > 0)
            {
                if (!select_id)
                {
                    break_dist <- 0.4
                    extra_cols <- setdiff(colnames(pdt),
                                          c("parameter", "np", "varying", "value"))
                    if (length(extra_cols) > 0)
                    {
                        cast_formula <-
                            stats::as.formula(paste(paste(c("np", extra_cols), collapse = " + "), "parameter", sep = "~"))
                    } else
                    {
                        cast_formula <- as.formula("np~parameter")
                    }
                    wpdt <-
                      data.table::data.table(data.table::dcast(pdt[varying == TRUE & distribution == "posterior"], cast_formula, value.var = "value"))
                    wpdt[, np := NULL]
                    if (length(extra_cols) > 0)
                    {
                        wpdt[, paste(extra_cols) := NULL]
                    }
                    cp <- GGally::ggcorr(wpdt)
                    plots[["correlations"]] <- cp
                    pp <- GGally::ggpairs(wpdt)
                    plots[["pairs"]] <- pp
                }

                density_data <- pdt[varying == TRUE]
                if (black_prior) {
                    density_data <- density_data[distribution == "posterior"]
                }
                dp <- ggplot()
                dp <- dp + facet_wrap(~ parameter, scales = "free",
                                      labeller = label_parsed)
                dp <- dp + do.call(paste0("geom_", densities), c(list(mapping = do.call(aes_string, aesthetic), data = density_data, position = "identity"),
                                                                 density_args))
                if (black_prior) {
                    dp <- dp + geom_line(data = pdt[varying == TRUE & distribution == "prior"], mapping = aes(x = value), stat = "density", color = "black", adjust = 2)
                }
                if (!missing(brewer.palette))
                {
                    dp <- dp + scale_color_brewer(palette = brewer.palette)
                    dp <- dp + scale_fill_brewer(palette = brewer.palette)
                }
                dp <- dp + ylab("density")
                dp <- dp + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                 legend.position = "top")
                if (select_id)
                {
                    dp <- dp + geom_vline(data = pdt[np %in% select[["np"]]],
                                          aes(xintercept = value))
                }
                plots[["densities"]] <- dp

                aesthetic <- list(x = "np", y = "value")

                if (!missing(extra.aes))
                {
                    aesthetic <- c(aesthetic, extra.aes)
                }

                tp <- ggplot(mapping = do.call(aes_string, aesthetic))
                tp <- tp + geom_line(data = pdt[varying == TRUE & distribution == "posterior"])
                tp <- tp + facet_wrap(~ parameter, scales = "free_y",
                                      labeller = label_parsed)
                tp <- tp + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                 legend.position = "top")
                if (select_id)
                {
                    tp <- tp + geom_vline(xintercept = select[["np"]])
                }
                if (!missing(brewer.palette))
                {
                    tp <- tp + scale_color_brewer(palette = brewer.palette)
                    tp <- tp + scale_fill_brewer(palette = brewer.palette)
                }
                plots[["traces"]] <- tp
            }
        }
    }

    if ("logeval" %in% type)
    {
        ldt <- NULL
        if (length(vars[["logeval"]]) > 0) {
            if (verbose) message(date(), " Getting logevals")
            samples <- clean_data(x, "x", verbose=verbose, vars=vars[["logeval"]])
        }

        if (verbose) message(date(), " Logeval plots")
        for (ll in vars[["logeval"]])
        {
            values <- samples[[ll]]
            if ("np" %in% colnames(samples[[ll]]))
            {
                values <- values[np >= burn]
            }

            if (!("data.frame" %in% class(values)))
            {
                values <- data.table::data.table(np = 0, value = values)
            }

            if (!("np" %in% colnames(values)))
            {
                data.table::setnames(values, "time", "np")
            }

            values[, density := ll]
            if (is.null(ldt))
            {
                ldt <- values
            } else
            {
                ldt <- rbind(ldt, values)
            }
        }
        trace_aesthetic <- list(x = "np", y = "value")
        density_aesthetic <- list(x = "value")

        if (!is.null(ldt) && nrow(ldt) > 0)
        {
            ret_data <- c(ret_data, list(likelihoods = ldt))

            likelihood_dummy_plot <-
                data.frame(expand.grid(type = c("Density", "Trace"),
                                       density = vars[["logeval"]]))
            lp <- ggplot(likelihood_dummy_plot)
            lp <- lp + geom_line(data = data.frame(ldt, type = "Trace"),
                                 mapping = do.call(aes_string, trace_aesthetic))
            lp <- lp +
                geom_histogram(data = data.frame(ldt, type = "Density"),
                               mapping = do.call(aes_string, density_aesthetic))
            lp <- lp + facet_wrap(density ~ type, scales = "free")
            lp <- lp + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                             legend.position = "top")
            if (select_id)
            {
                lp <- lp + geom_vline(data = data.frame(type = "Trace"),
                                      xintercept = ldt[np %in% select[["np"]], value])
            }
            plots[["logevals"]] <- lp
        }
    }

    ## plot first plot if requested
    if (plot)
    {
        print(plots[[1]])
    }
    plots[["data"]] <- ret_data

    ## return all plots invisibly
    invisible(plots)
}

##' Plot routing for \code{libbi} objects
##'
##' @param x \code{libbi} object
##' @param ... parameters to \code{\link{plot_libbi}}
##' @return a list of plots plot (see \code{\link{plot_libbi}})
##' @export
plot.libbi <- function(x, ...)
{
    plot_libbi(x, ...)
}
