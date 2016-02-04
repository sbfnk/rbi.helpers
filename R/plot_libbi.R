##' Plot results from libbi
##'
##' @param read Monte-Carlo samples, either a \code{libbi} object or a list of data frames, as returned by \code{bi_read}
##' @param prior optional; Prior samples, either a \code{libbi} object or a list of data frames, as returned by \code{bi_read}
##' @param model model file or a \code{bi_model} object (if \code{read} is not a \code{libbi} object)
##' @param states states to plot (if not given, all states will be plotted; if empty vector is passed, no states are plotted)
##' @param params parameters to plot (if not given, all states will be plotted; if empty vector is passed, no parameters are plotted)
##' @param noises noises to plot (if not given, all noises will be plotted; if empty vector is passed, no noises are plotted)
##' @param quantile.span if plots are produced, which quantile to use for confidence intervals
##' @param date.origin date of origin (if dates are to be calculated)
##' @param date.unit unit of date (if desired, otherwise the time dimension will be used)
##' @param time.dim time dimension ("nr" by default)
##' @param data data (with a "time" and "value" column)
##' @param id one or more run ids to plot
##' @param extra.aes extra aesthetics (for ggplot)
##' @param all.times whether to plot all times (not only ones with data)
##' @param hline horizontal marker lines, named vector in format (state = value)
##' @param burn How many runs to burn
##' @param thin How many runs to thin per run kept
##' @param steps whether to plot lines as stepped lines
##' @param select list of selection criteria
##' @param shift list of dimensions to be shifted, and by how much
##' @param data.colour colour for plotting the data
##' @param base.alpha base alpha value for credible intervals
##' @param trend how the trend should be characterised (e.g., mean, median)
##' @param densities density geometry (e.g., "histogram")
##' @param density_args list of arguments to pass to density geometry
##' @param limit.to.data whether to limit the time axis to times in the data
##' @param brewer.palette optional; brewer color palette
##' @param ... options for geom_step / geom_line
##' @return list of results
##' @import ggplot2 scales reshape2
##' @importFrom lubridate wday %m+% years
##' @export
##' @author Sebastian Funk
plot_libbi <- function(read, prior, model, states, params, noises,
                       quantile.span = c(0.5, 0.95),
                       date.origin, date.unit, time.dim = "nr",
                       data, id, extra.aes,
                       all.times = FALSE, hline,
                       burn, thin, steps = FALSE, select,
                       shift, data.colour = "red", base.alpha = 0.5,
                       trend = "median", densities = "density",
                       density_args = NULL, limit.to.data = FALSE,
                       brewer.palette, ...)
{
    use_dates <- FALSE
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

    ret_data <- list()
    ## copy data table

    if ("libbi" %in% class(read))
    {
        if (!read$run_flag)
        {
            stop("The model should be run first")
        }
        res <- bi_read(read)
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
    res <- lapply(res, function(x) { if (is.data.frame(x)) { data.table(x) } else {x} })
    res <- lapply(res, copy)

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
    }
    res_prior <- lapply(res_prior, function(x) { if (is.data.frame(x)) { data.table(x) } else {x} })
    res_prior <- lapply(res_prior, copy)

    if (missing(model))
    {
        if ("libbi" %in% class(read))
        {
            model <- read$model
        } else
        {
            stop("If 'read' is not a libbi object, 'model' must be given.")
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
    }

    if (!missing(burn) || !missing(thin))
    {
        for (var in names(res))
        {
            if ("np" %in% colnames(res[[var]]))
            {
                if (!missing(burn))
                {
                    res[[var]] <- res[[var]][np >= burn]
                }
                if (!missing(thin))
                {
                    res[[var]] <- res[[var]][np %% (thin + 1) == 0]
                }
                recount_np <- data.table(np = unique(res[[var]][, np]))
                recount_np[, new_np := seq_len(nrow(recount_np)) - 1]
                res[[var]] <- merge(res[[var]], recount_np, by = "np")
                res[[var]][, np := NULL]
                setnames(res[[var]], "new_np", "np")
            }
        }
    }

    sdt <- data.table(state = character(0))
    if (use_dates)
    {
        sdt[, time := as.Date(character(0))]
        sdt[, time_next := as.Date(character(0))]
    } else
    {
        sdt[, time := numeric(0)]
        sdt[, time_next := numeric(0)]
    }
    sdt[, value := numeric(0)]
    if (missing(id))
    {
        for (i in seq_along(quantile.span))
        {
            sdt[, paste("min", i, sep = ".") := numeric(0)]
            sdt[, paste("max", i, sep = ".") := numeric(0)]
        }
    } else
    {
        sdt[, np := numeric(0)]
    }
    if (!missing(extra.aes))
    {
        for (extra in unique(unname(extra.aes)))
        {
            sdt[, paste(extra) := character(0)]
        }
    }

    p <- NULL
    if (missing(states))
    {
        states <- c(model$get_vars("state"), model$get_vars("obs"))
    }

    states <- intersect(names(res), states)

    if (length(states) > 0)
    {
        if (!missing(data))
        {
            if (length(setdiff(c("time", "value"), colnames(data))) > 0) {
                stop("'data' does not have a 'time' and 'value' column.")
            } else {
                dataset <- data.table(data)
            }
        }

        for (state in states)
        {
            if (state %in% names(res))
            {
                values <- res[[state]]

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

                summarise_columns <- c("nr", "np", "time", "time_next")

                if (!missing(extra.aes))
                {
                    summarise_columns <- c(summarise_columns, unique(unname(extra.aes)))
                }
                sum.by <- intersect(summarise_columns, colnames(values))
                values <- values[, list(value = sum(value)), by = sum.by]

                state.by <- intersect(setdiff(summarise_columns, "np"),
                                      colnames(values))
                state.wo <- setdiff(setdiff(summarise_columns, "np"),
                                    colnames(values))

                if (missing(id))
                {
                    temp_values <-
                        values[, list(value = do.call(trend, list(value, na.rm = TRUE))),
                               by = state.by]
                    for (i in seq_along(quantile.span))
                    {
                        quantiles <-
                            values[, list(max = quantile(value, 0.5 + quantile.span[i] / 2, na.rm = TRUE),
                                          min = quantile(value, 0.5 - quantile.span[i] / 2, na.rm = TRUE)),
                                   by = state.by]
                        setnames(quantiles, c("min", "max"), paste(c("min",  "max"),  i,  sep = "."))
                        temp_values <- merge(temp_values, quantiles, by = state.by)
                    }
                    values <- temp_values
                    if (steps)
                    {
                        max.nr <- values[, max(nr)]
                        for (i in seq_along(quantile.span))
                        {
                            values[nr == max.nr, paste("min", i, sep = ".") := NA]
                            values[nr == max.nr, paste("max", i, sep = ".") := NA]
                        }
                    }
                } else
                {
                    if (!("all" %in% id))
                    {
                        values <- values[np %in% id]
                    }
                }

                for (wo in state.wo)
                {
                    values[, paste(wo) := "n/a"]
                }
                values[, paste(time.dim) := NULL]
                sdt <- rbind(sdt,
                             data.table(state = rep(state, nrow(values)),
                                        values))
            } else if (state != "dummy")
            {
                warning(paste("State", state, "does not exist"))
            }
        }
        sdt[, state := factor(state, levels = unique(state))]
        ret_data <- c(ret_data, list(states = sdt))

        if (!missing(data) && nrow(dataset) > 0)
        {
            dataset <- dataset[state %in% states]
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
                    for (data_state in unique(dataset[, state]))
                    {
                        if (limit.to.data)
                        {
                            ## for all states, only retain times in data
                            sdt <- sdt[time %in% dataset[state == data_state, time]]
                        } else
                        {
                            ## for states in dataset, only retain times in data
                            sdt <- sdt[(state != data_state) |
                                       (time %in% dataset[state == data_state, time])]
                        }
                    }
                }
                for (i in seq_along(quantile.span))
                {
                    dataset[, paste("min", i, sep = ".") := 0]
                    dataset[, paste("max", i, sep = ".") := 0]
                }
                ret_data <- c(ret_data, list(data = dataset))
            }
        }

        aesthetic <- list(x = "time", y = "value")
        if (!missing(id) && ("all" %in% id || length(id) > 1))
        {
            aesthetic <- c(aesthetic, list(color = "factor(np)"))
        }
        if (!missing(extra.aes))
        {
            aesthetic <- c(aesthetic, extra.aes)
        }

        if (nrow(sdt) > 0)
        {
            p <- ggplot(sdt, do.call(aes_string, aesthetic))

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
                                    ncol = round(sqrt(length(states))))
            }
            if (missing(id))
            {
                alpha <- base.alpha
                for (i in seq_along(quantile.span))
                {
                    str <- as.list(paste(c("max", "min"), i, sep = "."))
                    names(str) <- c("ymax", "ymin")
                    p <- p + ribbon_func(do.call(aes_string, str), alpha = alpha)
                    alpha <- alpha / 2
                }
            }
            if ("color" %in% names(aesthetic))
            {
                p <- p + line_func(...)
            } else
            {
                p <- p + line_func(color = "black", ...)
            }
            p <- p + scale_y_continuous("", labels = comma)
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

    pdt <- data.table(distribution = character(0),
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
        params <- model$get_vars("param")
    }

    params <- intersect(names(res), params)

    if (length(params) > 0)
    {
        for (param in params)
        {
            param_values <- list()
            param_values[["posterior"]] <- res[[param]]
            if (!is.null(res_prior) && param %in% names(res_prior))
            {
                param_values[["prior"]] <- res_prior[[param]]
            }

            for (dist in names(param_values))
            {
                values <- param_values[[dist]]
                if (!("data.frame" %in% class(values)))
                {
                    values <- data.table(np = 0, value = values)
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
                             data.table(distribution = dist,
                                        parameter = rep(sub("^p_", "", param),
                                                        nrow(values)),
                                        values))
            }
        }
        if (nrow(pdt) > 0)
        {
            pdt[, parameter := factor(parameter, levels = unique(parameter))]

            by.varying <- c("parameter", "distribution")
            if (!missing(extra.aes))
            {
                by.varying <- c(by.varying, unique(unname(extra.aes)))
            }
            pdt[, varying := (length(unique(value)) > 1), by = by.varying]

            ret_data <- c(ret_data, list(params = pdt))

            aesthetic <- list(x = "value")
            if (!missing(extra.aes))
            {
                aesthetic <- c(aesthetic, extra.aes)
            }

            if (!is.null(res_prior))
            {
                if (!missing(extra.aes) && length(intersect(c("color", "fill"), names(extra.aes))) > 0)
                {
                    warning("'extra.aes' (with color or fill) and 'prior' given. Will ignore 'prior' as things will look a mess otherwise.")
                } else
                {
                    aesthetic <- c(aesthetic, list(color = "distribution", fill = "distribution"))
                }
            }

            if (nrow(pdt[varying == TRUE]) > 0)
            {
                if (missing(id))
                {
                    break_dist <- 0.4
                    extra_cols <- setdiff(colnames(pdt),
                                          c("parameter", "np", "varying", "value"))
                    if (length(extra_cols) > 0)
                    {
                        cast_formula <-
                            as.formula(paste(paste("np", extra_cols, sep = "+"),
                                             "parameter", sep = "~"))
                    } else
                    {
                        cast_formula <- as.formula("np~parameter")
                    }
                    wpdt <- data.table(dcast(pdt[varying == TRUE & distribution == "posterior"],
                                             cast_formula,
                                             value.var = "value"))
                    wpdt[, np := NULL]
                    if (length(extra_cols) > 0)
                    {
                        wpdt[, paste(extra_cols) := NULL]
                    }
                    correlations <- data.table(melt(cor(wpdt)))
                    correlations[, correlation := cut(value,
                                                      breaks = c(seq(-1, 1, by = break_dist)))]
                    ret_data <- c(ret_data, list(correlations = correlations))

                    color_palette <-
                        colorRampPalette(c("#3794bf", "#FFFFFF", "#df8640"))(2 / break_dist)

                    cp <- ggplot(correlations, aes(x = Var1, y = Var2,
                                                   fill = correlation))
                    cp <- cp + geom_tile()
                    cp <- cp + scale_fill_manual("Correlation", values = color_palette,
                                                 limits = levels(correlations[, correlation]))
                    cp <- cp + scale_x_discrete("")
                    cp <- cp + scale_y_discrete("")
                } else {
                    cp <- NULL
                }

                dp <- ggplot(pdt[varying == TRUE], do.call(aes_string, aesthetic))
                dp <- dp + facet_wrap(~ parameter, scales = "free")
                dp <- dp + do.call(paste0("geom_", densities),
                                   c(list(alpha = 0.5), density_args))
                if (!missing(brewer.palette))
                {
                    dp <- dp + scale_color_brewer(palette = brewer.palette)
                    dp <- dp + scale_fill_brewer(palette = brewer.palette)
                }
                dp <- dp + scale_y_continuous("Density")
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

                tp <- ggplot(pdt[varying == TRUE & distribution == "posterior"],
                             do.call(aes_string, aesthetic))
                tp <- tp + geom_line()
                tp <- tp + facet_wrap(~ parameter, scales = "free_y")
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

    ndt <- data.table(noise = character(0))
    if (use_dates)
    {
        ndt[, time := as.Date(character(0))]
        ndt[, time_next := as.Date(character(0))]
    } else
    {
        ndt[, time := numeric(0)]
        ndt[, time_next := numeric(0)]
    }
    ndt[, value := numeric(0)]
    if (missing(id))
    {
        for (i in seq_along(quantile.span))
        {
            ndt[, paste("min", i, sep = ".") := numeric(0)]
            ndt[, paste("max", i, sep = ".") := numeric(0)]
        }
    }
    if (!missing(extra.aes))
    {
        for (extra in unique(unname(extra.aes)))
        {
            ndt[, paste(extra) := character(0)]
        }
    }

    np <- NULL

    if (missing(noises))
    {
        noises <- model$get_vars("noise")
    }

    noises <- intersect(names(res), noises)

    if (length(noises) > 0)
    {
        for (noise in noises)
        {
            values <- res[[noise]]
            values[!is.finite(value), value := 0]

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

            by.sum <- c("nr", "time", "time_next")
            if (!missing(extra.aes))
            {
                by.sum <- c(by.sum, unique(unname(extra.aes)))
            }
            noise.by <- intersect(by.sum, colnames(values))
            noise.wo <- setdiff(by.sum, colnames(values))

            if (missing(id))
            {
                temp_values <-
                    values[, list(value = do.call(trend, list(value, na.rm = TRUE))),
                           by = noise.by]
                for (i in seq_along(quantile.span))
                {
                    quantiles <-
                        values[, list(max = quantile(value, 0.5 + quantile.span[i] / 2, na.rm = TRUE),
                                      min = quantile(value, 0.5 - quantile.span[i] / 2, na.rm = TRUE)),
                               by = noise.by]
                    setnames(quantiles, c("min", "max"), paste(c("min",  "max"),  i,  sep = "."))
                    temp_values <- merge(temp_values, quantiles, by = noise.by)
                }
                values <- temp_values
            } else
            {
                values <- values[np %in% id]
                values <- values[, list(value = sum(value)), by = noise.by]
            }
            for (wo in noise.wo)
            {
                values[, paste(wo) := "n/a"]
            }
            values[, paste(time.dim) := NULL]
            ndt <- rbind(ndt,
                         data.table(noise = rep(noise, nrow(values)),
                                    values))

            ndt[, noise := factor(noise, levels = unique(noise))]

            if (!missing(data) && nrow(dataset) > 0 && !all.times)
            {
                ndt <- ndt[(time >= min(dataset[, time])) & (time <= max(dataset[, time]))]
            }

            ret_data <- c(ret_data, list(noises = ndt))

            aesthetic <- list(x = "time", y = "value")
            if (!missing(extra.aes))
            {
                aesthetic <- c(aesthetic, extra.aes)
            }

            if (nrow(ndt) > 0)
            {
                np <- ggplot(ndt, do.call(aes_string, aesthetic))
                np <- np + facet_wrap(~ noise, scales = "free_y")
                if (missing(id))
                {
                    alpha <- base.alpha
                    for (i in seq_along(quantile.span))
                    {
                        str <- as.list(paste(c("max", "min"), i, sep = "."))
                        names(str) <- c("ymin", "ymax")
                        np <- np + geom_ribbon(do.call(aes_string, str), alpha = alpha)
                        alpha <- alpha / 2
                    }
                }
                np <- np + geom_line()
                np <- np + scale_y_continuous("", labels = comma)
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
    }

    lp <- NULL
    likelihoods <- intersect(names(res), c("loglikelihood", "logprior"))
    if (length(likelihoods) > 0)
    {
        ldt <- NULL
        for (ll in likelihoods)
        {
            values <- res[[ll]]

            if (!("data.frame" %in% class(values)))
            {
                values <- data.table(np = 0, value = values)
            }

            if (!("np" %in% colnames(values)))
            {
                setnames(values, "nr", "np")
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
        extra_columns <- setdiff(colnames(ldt), c("np", "value", "density"))
        if (length(extra_columns) > 0)
        {
            temp_ldt <-
                ldt[, list(value = do.call(trend, list(value, na.rm = TRUE))),
                    by = c("np", "density")]
            for (i in seq_along(quantile.span))
            {
                quantiles <-
                    ldt[, list(max = quantile(value, 0.5 + quantile.span[i] / 2, na.rm = TRUE),
                               min = quantile(value, 0.5 - quantile.span[i] / 2, na.rm = TRUE)),
                        by = c("np", "density")]
                setnames(quantiles, c("min", "max"), paste(c("min",  "max"),  i,  sep = "."))
                temp_ldt <- merge(temp_ldt, quantiles, by = c("np", "density"))
            }
            ldt <- temp_ldt
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
            if (length(extra_columns) > 0)
            {
                alpha <- base.alpha
                for (i in seq_along(quantile.span))
                {
                    str <- as.list(paste(c("max", "min"), i, sep = "."))
                    names(str) <- c("ymin", "ymax")
                    lp <- lp + geom_ribbon(do.call(aes_string, str), alpha = alpha)
                    alpha <- alpha / 2
                }
            }
            lp <- lp +
                geom_histogram(data = data.frame(ldt, type = "Density"),
                               mapping = do.call(aes_string, density_aesthetic))
            lp <- lp + facet_wrap(density ~ type, scales = "free")
            lp <- lp + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                             legend.position = "top")
            if (!missing(id))
            {
                lp <- lp + geom_vline(data = data.frame(type = "Trace"),
                                      xintercept = id)
            }
        }
    }

    return(list(states = p,
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
##' @return plot
plot.libbi <- function(obj, ...)
{
    plot_libbi(obj, ...)
}
