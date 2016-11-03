##' Plot results from libbi
##'
##' Plots state trajectories (unless plot = FALSE) and invisibly returns a list of state trajectories and other plots.
##' @param read Monte-Carlo samples, either a \code{libbi} object or a list of data frames, as returned by \code{bi_read}
##' @param model model file or a \code{bi_model} object (if \code{read} is not a \code{libbi} object)
##' @param prior optional; Prior samples, either a \code{libbi} object or a list of data frames, as returned by \code{bi_read}
##' @param states states to plot (if not given, all states will be plotted; if empty vector is passed, no states are plotted)
##' @param params parameters to plot (if not given, all states will be plotted; if empty vector is passed, no parameters are plotted)
##' @param noises noises to plot (if not given, all noises will be plotted; if empty vector is passed, no noises are plotted)
##' @param quantiles if plots are produced, which quantile to use for confidence intervals (NULL for no confidence intervals)
##' @param date.origin date of origin (if dates are to be calculated)
##' @param date.unit unit of date (if desired, otherwise the time dimension will be used)
##' @param time.dim time dimension ("time" by default)
##' @param data data (with a "time" and "value" column)
##' @param id one or more run ids to plot
##' @param extra.aes extra aesthetics (for ggplot)
##' @param all.times whether to plot all times (not only ones with data)
##' @param hline horizontal marker lines, named vector in format (state = value)
##' @param burn How many iterations to burn
##' @param steps whether to plot lines as stepped lines
##' @param select list of selection criteria
##' @param shift list of dimensions to be shifted, and by how much
##' @param data.colour colour for plotting the data
##' @param base.alpha base alpha value for credible intervals
##' @param trend how the trend should be characterised (e.g., mean, median, or NULL for no trend line)
##' @param densities density geometry (e.g., "histogram" (default) or "density")
##' @param density_args list of arguments to pass to density geometry
##' @param limit.to.data whether to limit the time axis to times in the data
##' @param labels facet labels, in case they are to be rewritten, to be parsed using \code{label_parsed}; should be given as named character vector of (parameter = 'label') pairs
##' @param brewer.palette optional; brewer color palette
##' @param plot set to FALSE to suppress plot of trajectories
##' @param ... options for geom_step / geom_line / geom_point / etc.
##' @return a list of plots
##' @import ggplot2 scales
##' @importFrom lubridate wday %m+% years
##' @importFrom rbi bi_read
##' @importFrom reshape2 dcast
##' @importFrom GGally ggcorr
##' @export
##' @author Sebastian Funk
plot_libbi <- function(read, model, prior, states, params, noises,
                       quantiles = c(0.5, 0.95),
                       date.origin, date.unit, time.dim = "time",
                       data, id, extra.aes,
                       all.times = FALSE, hline,
                       burn, steps = FALSE, select,
                       shift, data.colour = "red", base.alpha = 0.5,
                       trend = "median", densities = "histogram",
                       density_args = list(), limit.to.data = FALSE,
                       labels, brewer.palette, plot = TRUE, ...)
{
    use_dates <- FALSE
    summarise_columns <- c("np", "time", "time_next")

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

    ret_data <- list()
    ## copy data table

    if ("libbi" %in% class(read))
    {
        if (!read$run_flag)
        {
            stop("The model should be run first")
        }
        res <- rbi::bi_read(read)
    } else if (is.data.frame(read))
    {
        res <- list(dummy = read)
    } else if (is.list(read))
    {
        res <- read
    } else
    {
        stop("'read' must be a 'libbi' object or a list of data frames or a data frame.")
    }
    res <- lapply(res, function(x) { if (is.data.frame(x)) { data.table::data.table(x) } else {x} })

    res_prior <- NULL
    if (!missing(prior))
    {
        if ("libbi" %in% class(prior))
        {
            if (!prior$run_flag)
            {
                stop("The model should be run first")
            }
            res_prior <- bi_read(prior)
        } else if (is.data.frame(prior))
        {
            res_prior <- list(dummy = prior)
        } else if (is.list(prior))
        {
            res_prior <- prior
        } else
        {
            stop("'prior' must be a 'libbi' object or a list of data frames or a data frame.")
        }
        res_prior <- lapply(res_prior, function(x) { if (is.data.frame(x)) { data.table::data.table(x) } else {x} })
    }

    if (missing(model))
    {
        if ("libbi" %in% class(read))
        {
            model <- read$model
        }
    } else
    {
        if ("libbi" %in% class(read))
        {
            warning("'model' overwrites the model given in 'read'.x")
        }
        if (is.character(model)) {
            model <- bi_model(model)
        } else if (!("bi_model" %in% class(model))) {
            stop("'model' must be either a 'bi_model' object or a path to a valid model file in LibBi's syntax")
        }
    }

    if (steps)
    {
        ribbon_func <- function(mapping, ...)
        {
            if (missing(mapping)) mapping <- list()
            mapping <- modifyList(mapping, aes_string(xmin = "time"))
            mapping <- modifyList(mapping, aes_string(xmax = "time_next"))
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

    sdt <- NULL

    p <- NULL
    if (missing(states))
    {
        if (missing(model))
        {
            states <- c()
            given_states <- c()
        } else
        {
            states <- c(model$get_vars("state"), model$get_vars("obs"))
            given_states <- c()
        }
    } else
    {
        given_states <- states
    }

    states <- intersect(names(res), states)

    missing_states <- setdiff(given_states, states)
    if (length(missing_states) > 0)
    {
        warning("State(s) ", missing_states, " not found in data.")
    }

    if (length(states) > 0)
    {
        if (!missing(data))
        {
            if (length(setdiff(c("time", "value"), colnames(data))) > 0) {
                stop("'data' does not have a 'time' and 'value' column.")
            } else {
                dataset <- data.table::data.table(data)
            }
        }

        for (state in states)
        {
          if (state %in% names(res) && nrow(res[[state]]) > 0)
            {
                values <- res[[state]]
                if ("np" %in% colnames(res[[state]]))
                {
                    values <- values[np >= burn]
                    if (nrow(values) == 0) {
                        stop("Nothing left after burn-in")
                    }
                }

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

                if (!missing(select))
                {
                    for (var_name in names(select))
                    {
                        if (var_name %in% colnames(values))
                        {
                            values <- values[get(var_name) %in% select[[var_name]]]
                            if (class(values[, get(var_name)]) == "factor")
                            {
                                values[, paste(var_name) := factor(get(var_name))]
                            }
                        }
                    }
                }

                if (!missing(shift))
                {
                    for (var_name in names(shift))
                    {
                        if (var_name %in% colnames(values))
                        {
                            values <- values[get(var_name) >= shift[[var_name]]]
                            values[, paste(var_name) := get(var_name) - shift[[var_name]]]
                        }
                    }
                }

                sum.by <- intersect(summarise_columns, colnames(values))
                values <- values[, list(value = sum(value)), by = sum.by]

                state.wo <- setdiff(setdiff(summarise_columns, "np"),
                                    colnames(values))
                for (wo in state.wo)
                {
                    values[, paste(wo) := "n/a"]
                }

                new_states <- data.table::data.table(state = rep(state, nrow(values)), values)
                if (is.null(sdt))
                {
                    sdt <- new_states
                } else
                {
                    sdt <- rbind(sdt, new_states)
                }
            } else if (state != "dummy")
            {
                warning(paste("State", state, "does not exist"))
            }
        }
        ## factorise columns
        sdt <- factorise_columns(sdt, labels)

        if (!missing(data) && nrow(dataset) > 0)
        {
            if (nrow(dataset) > 0)
            {
                if (!missing(select))
                {
                    for (var_name in names(select))
                    {
                        if (var_name %in% unique(sdt[, state]))
                        {
                            dataset <- dataset[get(var_name) %in% select[[var_name]]]
                            if (class(values[, get(var_name)]) == "factor")
                            {
                                dataset[, paste(var_name) := factor(get(var_name))]
                            }
                        }
                    }
                }

                if (!all.times && nrow(sdt) > 0)
                {
                    if (limit.to.data)
                    {
                        ## for all states, only retain times in data
                        sdt <- sdt[time %in% dataset[, time]]
                    } else
                    {
                        for (data_state in unique(dataset[, state]))
                        {
                            ## for states in dataset, only retain times in data
                            sdt <- sdt[(state != data_state) |
                                       (time %in% dataset[state == data_state, time])]
                        }
                    }
                }
                for (i in seq_along(quantiles))
                {
                    dataset[, paste("min", i, sep = ".") := 0]
                    dataset[, paste("max", i, sep = ".") := 0]
                }

                for (missing_column in setdiff(colnames(sdt), colnames(dataset)))
                {
                    dataset[, paste(missing_column) := "n/a"]
                }
                ret_data <- c(ret_data, list(data = dataset))
            }
        }

        states_n <- sdt[, list(single = (.N == 1)), by = state]
        sdt <- merge(sdt, states_n, by = "state", all.x = TRUE)

        if (!missing(id) && !("all" %in% id))
        {
            sdt <- sdt[np %in% id]
        }

        aggregate_values <- NULL
        state.by <- c("state", "single", intersect(setdiff(summarise_columns, "np"), colnames(sdt)))
        if (is.null(trend))
        {
            id_alpha = 1
        } else
        {
            aggregate_values <- sdt[, list(value = do.call(trend, list(value, na.rm = TRUE))),
                                       by = state.by]
            id_alpha = 0.35
        }

        if (!is.null(quantiles))
        {
            for (i in seq_along(quantiles))
            {
                quantile_values <-
                    sdt[, list(max = quantile(value, 0.5 + quantiles[i] / 2, na.rm = TRUE),
                                  min = quantile(value, 0.5 - quantiles[i] / 2, na.rm = TRUE)),
                        by = state.by]
                setnames(quantile_values, c("min", "max"), paste(c("min",  "max"),  i,  sep = "."))
                if (is.null(aggregate_values))
                {
                    aggregate_values <- quantile_values
                } else
                {
                    aggregate_values <- merge(aggregate_values, quantile_values, by = state.by)
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
        if (is.null(aggregate_values)) aggregate_values <- sdt

        ret_data <- c(ret_data, list(states = aggregate_values[, !"single", with = FALSE]))

        aesthetic <- list(x = "time", y = "value")
        if (!missing(extra.aes))
        {
            aesthetic <- c(aesthetic, extra.aes)
        }

        if (nrow(aggregate_values) > 0)
        {
            if (!missing(extra.aes) && "color" %in% names(extra.aes) &&
                !missing(id))
            {
                sdt <- sdt[, color_np := paste(get(extra.aes["color"]), get("np"), sep = "_")]
            }
            p <- ggplot(mapping = do.call(aes_string, aesthetic))

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
                for (hline_state_id in named)
                {
                    hline_data <- data.frame(state = names(hline)[hline_state_id],
                                             yintercept = hline[hline_state_id])
                    p <- p + geom_hline(data = hline_data,
                                        aes(yintercept = yintercept), color = "black")
                }
                for (hline_state_id in unnamed)
                {
                    hline_data <- data.frame(yintercept = hline[hline_state_id])
                    p <- p + geom_hline(data = hline_data,
                                        aes(yintercept = yintercept), color = "black")
                }
            }

            if (length(states) > 1)
            {
                p <- p + facet_wrap(~ state, scales = "free_y",
                                    ncol = round(sqrt(length(states))),
                                    labeller = label_parsed)
            }
            if (!is.null(quantiles))
            {
                alpha <- base.alpha
                for (i in seq_along(quantiles))
                {
                    str <- as.list(paste(c("max", "min"), i, sep = "."))
                    names(str) <- c("ymax", "ymin")
                    p <- p + ribbon_func(data = aggregate_values, do.call(aes_string, str), alpha = alpha)
                    alpha <- alpha / 2
                    if (nrow(sdt[single == TRUE]) > 0)
                    {
                        p <- p + geom_errorbar(data = aggregate_values[single == TRUE], do.call(aes_string, str), ...)
                    }
                }
            }
            if ("color" %in% names(aesthetic))
            {
                if (!is.null(trend))
                {
                    p <- p + line_func(data = aggregate_values, ...)
                    if (nrow(aggregate_values[single == TRUE]) > 0)
                    {
                        p <- p + geom_point(data = aggregate_values[single == TRUE], shape = 4, ...)
                    }
                }
                if (!missing(id))
                {
                    p <- p + line_func(data = sdt[single == FALSE], mapping = aes(group = color_np), alpha = id_alpha, ...)
                    if (nrow(sdt[single == TRUE]) > 0)
                    {
                        p <- p + geom_point(data = aggregate_values[single == TRUE], aes(group = factor(np)), shape = 4, alpha = id_alpha, ...)
                    }
                }
            } else
            {
                if (!is.null(trend))
                {
                    p <- p + line_func(data = aggregate_values, ...)
                    p <- p + geom_point(data = aggregate_values[single == TRUE], shape = 4, ...)
                }
                if (!missing(id))
                {
                    p <- p + line_func(data = sdt[single == FALSE], aes(group = factor(np)), alpha = id_alpha, ...)
                    p <- p + geom_point(data = sdt[single == TRUE], aes(group = factor(np)), alpha = id_alpha, shape = 4, ...)
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
        }
    }

    pdt <- data.table::data.table(distribution = character(0),
                                  parameter = character(0), np = integer(0),
                                  value = numeric(0))
    if (!missing(extra.aes))
    {
        for (extra in unique(unname(extra.aes)))
        {
            pdt[, paste(extra) := character(0)]
        }
    }

    if (missing(params))
    {
        if (missing(model))
        {
            params <- c()
        } else
        {
            params <- model$get_vars("param")
        }
    }

    params <- intersect(names(res), params)

    if (length(params) > 0)
    {
        for (param in params)
        {
            param_values <- list()
            param_values[["posterior"]] <- res[[param]]
            if ("np" %in% colnames(res[[param]]))
            {
                param_values[["posterior"]] <-
                    param_values[["posterior"]][np >= burn]
            }

            if (!is.null(res_prior) && param %in% names(res_prior))
            {
                param_values[["prior"]] <- res_prior[[param]]
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
                pdt <- rbind(pdt,
                             data.table::data.table(distribution = dist,
                                                    parameter = rep(param, nrow(values)),
                                                    values))
            }
        }
        if (nrow(pdt) > 0)
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
            if (!is.null(res_prior))
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

            if (nrow(pdt[varying == TRUE & distribution == "posterior"]) > 0)
            {
                if (missing(id))
                {
                    break_dist <- 0.4
                    extra_cols <- setdiff(colnames(pdt),
                                          c("parameter", "np", "varying", "value"))
                    if (length(extra_cols) > 0)
                    {
                        cast_formula <-
                            as.formula(paste(paste(c("np", extra_cols), collapse = " + "),
                                             "parameter", sep = "~"))
                    } else
                    {
                        cast_formula <- as.formula("np~parameter")
                    }
                    wpdt <-
                      data.table::data.table(reshape2::dcast(pdt[varying == TRUE & distribution == "posterior"],
                                                             cast_formula,
                                                             value.var = "value"))
                    wpdt[, np := NULL]
                    if (length(extra_cols) > 0)
                    {
                        wpdt[, paste(extra_cols) := NULL]
                    }
                    cp <- GGally::ggcorr(wpdt)
                } else {
                    cp <- NULL
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
                if (!missing(id))
                {
                    dp <- dp + geom_vline(data = pdt[np %in% id],
                                          aes(xintercept = value))
                }

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
                if (!missing(id))
                {
                    tp <- tp + geom_vline(xintercept = id)
                }
                if (!missing(brewer.palette))
                {
                    tp <- tp + scale_color_brewer(palette = brewer.palette)
                    tp <- tp + scale_fill_brewer(palette = brewer.palette)
                }
            } else
            {
                dp <- NULL
                tp <- NULL
                cp <- NULL
            }
        } else
        {
            dp <- NULL
            tp <- NULL
            cp <- NULL
        }
    } else
    {
        dp <- NULL
        tp <- NULL
        cp <- NULL
    }

    np <- NULL
    ndt <- NULL
    if (missing(noises))
    {
        if (missing(model))
        {
            noises <- c()
        } else
        {
            noises <- model$get_vars("noise")
        }
    }

    noises <- intersect(names(res), noises)

    if (length(noises) > 0)
    {
        for (noise in noises)
        {
            if (noise %in% names(res))
            {
              values <- res[[noise]]
                if ("np" %in% colnames(res[[noise]]))
                {
                    values <- values[np >= burn]
                }

              ## values[!is.finite(value), value := 0]

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

              sum.by <- intersect(summarise_columns, colnames(values))
              values <- values[, list(value = sum(value)), by = sum.by]

              state.wo <- setdiff(setdiff(summarise_columns, "np"),
                                  colnames(values))
              for (wo in state.wo)
              {
                values[, paste(wo) := "n/a"]
              }

              new_noises <- data.table::data.table(noise = rep(noise, nrow(values)), values)
              if (is.null(ndt))
              {
                ndt <- new_noises
              } else
              {
                ndt <- rbind(ndt, new_noises)
              }
            }
        }

        ndt <- factorise_columns(ndt, labels)

        noises_n <- ndt[, list(single = (.N == 1)), by = noise]
        ndt <- merge(ndt, noises_n, by = "noise", all.x = TRUE)

        aggregate_noises <- NULL
        noise.by <- c("noise", "single", intersect(setdiff(summarise_columns, "np"), colnames(ndt)))
        if (!is.null(trend))
        {
            aggregate_noises <-
                ndt[, list(value = do.call(trend,
                                           list(value, na.rm = TRUE))),
                    by = noise.by]
        }

        if (!is.null(quantiles))
        {
            for (i in seq_along(quantiles))
            {
                quantile_values <-
                    ndt[, list(max = quantile(value, 0.5 + quantiles[i] / 2, na.rm = TRUE),
                               min = quantile(value, 0.5 - quantiles[i] / 2, na.rm = TRUE)),
                        by = noise.by]
                setnames(quantile_values, c("min", "max"), paste(c("min",  "max"),  i,  sep = "."))
                if (is.null(aggregate_noises))
                {
                    aggregate_noises <- quantile_values
                } else
                {
                    aggregate_noises <- merge(aggregate_noises, quantile_values, by = noise.by)
                }
            }
            if (steps)
            {
                max.time <- aggregate_noises[, max(time)]
                for (i in seq_along(quantiles))
                {
                    aggregate_noises[time == max.time, paste("min", i, sep = ".") := NA]
                    aggregate_noises[time == max.time, paste("max", i, sep = ".") := NA]
                }
            }
        }

        if (!is.null(aggregate_noises))
        {
          ret_data <- c(ret_data, list(noises = aggregate_noises[, !"single", with = FALSE]))
        }

        if (!missing(id) && !("all" %in% id))
        {
            ndt <- ndt[np %in% id]
        }

        if (!missing(data) && nrow(dataset) > 0 && !all.times)
        {
            ndt <- ndt[(time >= min(dataset[, time])) & (time <= max(dataset[, time]))]
        }

        aesthetic <- list(x = "time", y = "value")
        if (!missing(extra.aes))
        {
            aesthetic <- c(aesthetic, extra.aes)
        }

        if (nrow(ndt) > 0)
        {
            if (!missing(extra.aes) && "color" %in% names(extra.aes) &&
                !missing(id))
            {
                ndt <- ndt[, color_np := paste(get(extra.aes["color"]), get("np"), sep = "_")]
            }
            np <- ggplot(mapping = do.call(aes_string, aesthetic))

            if (length(noises) > 1)
            {
                np <- np + facet_wrap(~ noise, scales = "free_y",
                                      ncol = round(sqrt(length(states))),
                                      labeller = label_parsed)
            }
            if (!is.null(quantiles))
            {
                alpha <- base.alpha
                for (i in seq_along(quantiles))
                {
                    str <- as.list(paste(c("max", "min"), i, sep = "."))
                    names(str) <- c("ymax", "ymin")
                    np <- np + ribbon_func(data = aggregate_noises, do.call(aes_string, str), alpha = alpha)
                    alpha <- alpha / 2
                    if (nrow(ndt[single == TRUE]) > 0)
                    {
                        np <- np + geom_errorbar(data = aggregate_noises[single == TRUE], do.call(aes_string, str), ...)
                    }
                }
            }
            if ("color" %in% names(aesthetic))
            {
                if (!is.null(trend))
                {
                    np <- np + line_func(data = aggregate_noises, ...)
                    if (nrow(aggregate_noises[single == TRUE]) > 0)
                    {
                        np <- np + geom_point(data = aggregate_noises[single == TRUE], shape = 4, ...)
                    }
                }
                if (!missing(id))
                {
                    np <- np + line_func(data = ndt[single == FALSE], mapping = aes(group = color_np), alpha = id_alpha, ...)
                    if (nrow(ndt[single == TRUE]) > 0)
                    {
                        np <- np + geom_point(data = ndt[single == TRUE], aes(group = factor(np)), shape = 4, alpha = id_alpha, ...)
                    }
                }
            } else
            {
                if (!is.null(trend))
                {
                    np <- np + line_func(data = aggregate_noises, ...)
                    np <- np + geom_point(data = aggregate_noises[single == TRUE], shape = 4, ...)
                }
                if (!missing(id))
                {
                    np <- np + line_func(data = ndt[single == FALSE], aes(group = factor(np)), alpha = id_alpha, ...)
                    np <- np + geom_point(data = ndt[single == TRUE], aes(group = factor(np)), alpha = id_alpha, shape = 4, ...)
                }
            }
            np <- np + scale_y_continuous(labels = comma) + ylab("")
            if (!missing(brewer.palette))
            {
                np <- np + scale_color_brewer(palette = brewer.palette)
                np <- np + scale_fill_brewer(palette = brewer.palette)
            }
            np <- np + expand_limits(y = 0)
            np <- np + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                              legend.position = "top")
            if (use_dates)
            {
                np <- np + scale_x_date("")
            }
        }
    }

    lp <- NULL
    likelihoods <- intersect(names(res), c("loglikelihood", "logprior"))
    if (length(likelihoods) > 0)
    {
        ldt <- NULL
        for (ll in likelihoods)
        {
            values <- res[[ll]]
            if ("np" %in% colnames(res[[ll]]))
            {
                values <- values[np >= burn]
            }

            if (!("data.frame" %in% class(values)))
            {
                values <- data.table::data.table(np = 0, value = values)
            }

            if (!("np" %in% colnames(values)))
            {
                setnames(values, "time", "np")
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

        ret_data <- c(ret_data, list(likelihoods = ldt))

        likelihood_dummy_plot <-
            data.frame(expand.grid(type = c("Density", "Trace"),
                                   density = c("loglikelihood", "logprior")))
        if (nrow(ldt) > 0)
        {
            lp <- ggplot(likelihood_dummy_plot)
            lp <- lp + geom_line(data = data.frame(ldt, type = "Trace"),
                                 mapping = do.call(aes_string, trace_aesthetic))
            lp <- lp +
                geom_histogram(data = data.frame(ldt, type = "Density"),
                               mapping = do.call(aes_string, density_aesthetic))
            lp <- lp + facet_wrap(density ~ type, scales = "free")
            lp <- lp + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                             legend.position = "top")
            if (!missing(id))
            {
              lp <- lp + geom_vline(data = data.frame(type = "Trace"),
                                    xintercept = ldt[np %in% id, value])
            }
        }
    }

    ## plot state trajectories unless told otherwise
    if (plot) print(p)

    ## return all plots invisibly
    invisible(list(states = p,
                   densities = dp,
                   traces = tp,
                   correlations = cp,
                   noises = np,
                   likelihoods = lp,
                   data = ret_data))
}

##' Plot routing for \code{libbi} objects
##'
##' @param obj \code{libbi} object
##' @param ... parameters to \code{\link{plot_libbi}}
##' @return a list of plots plot (see \code{\link{plot_libbi}})
##' @export
plot.libbi <- function(obj, ...)
{
    plot_libbi(obj, ...)
}
