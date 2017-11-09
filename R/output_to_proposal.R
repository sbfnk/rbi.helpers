#' @rdname output_to_proposal
#' @name output_to_proposal
#' @title Construct a proposal from run results
#' @description
#' This function takes the provided \code{\link{libbi}} which has been
#' run and returns a new model which has the proposal constructed from
#' the sample mean and standard deviation.
#' @param wrapper a \code{\link{libbi}} which has been run
#' @param scale a factor by which to scale all the standard deviations
#' @param correlations whether to take into account correlations
#' @param start whether this is the first attempt, in which case we'll use 1/10 of every bound, and 1 otherwise
#' @importFrom data.table setnames
#' @importFrom stats cov
#' @importFrom rbi get_block add_block insert_lines
#' @return the updated bi model
#' @keywords internal
output_to_proposal <- function(wrapper, scale, correlations = FALSE, start = FALSE) {

  if (!wrapper$run_flag) {
    stop("The model should be run first")
  }

  model <- wrapper$model
  ## get constant expressions
  const_lines <- grep("^[[:space:]]*const", model[], value = TRUE)
  for (const_line in const_lines) {
    line <-
      gsub(" ", "", sub("^[[:space:]]*const[[:space:]]*", "", const_line))
    assignment <- strsplit(line, "=")[[1]]
    tryCatch(
    {
      assign(assignment[1], eval(parse(text = assignment[2])))
    },
    error = function(cond)
    {
      assign(assignment[1], assignment[2])
    })
  }

  params <- list()
  param_block <- list()
  for (block in c("parameter", "initial"))
  {
    ## get parameters
    param_block[[block]] <- get_block(model, block)
    ## only go over variable parameters
    var_lines <- grep("~", param_block[[block]], value = TRUE)
    if (length(var_lines) > 0)
    {
      params[[block]] <- sub("[[:space:]]*~.*$", "", var_lines)
      ## remove any dimensions
      params[[block]] <- gsub("[[:space:]]*\\[[^]]*\\]", "", params[[block]])
    }
  }
  all_params <- unname(unlist(params))

  ## read parameters
  init_params <- ("with-transform-initial-to-param" %in% names(wrapper$options))
  res <- bi_read(wrapper, vars = all_params, init.to.param = init_params)

  if (correlations && length(all_params) > 1) {
    ## adapt to full covariance matrix
    if (!("__diff_" %in% var_names(model, "param", aux=TRUE))) {
      model <-
        rbi::insert_lines(model, "param __diff_ (has_output=0, has_input=0)",
                          after = 1)
    }

    l <- list()
    for (param in names(res)) {
      y <- copy(res[[param]])
      ## extract columns that are dimensions
      dim_colnames <- setdiff(colnames(y), c("np", "value"))
      unique_dims <- unique(y[dim_colnames])
      if (sum(dim(unique_dims)) > 0)
      {
        ## for parameters with dimensions, create a parameter for each
        ## possible dimension(s)
        a <- apply(unique_dims, 1, function(x) {
          merge(data.table(t(x)), y)
        })
        ## convert factor to integer
        for (dim_colname in dim_colnames) {
          if (is.factor(unique_dims[[dim_colname]])) {
            unique_dims[[dim_colname]] <-
              as.integer(unique_dims[[dim_colname]]) - 1
          }
        }
        ## create correct parameter names (including the dimensions)
        if (length(a))
          names(a) <- unname(apply(unique_dims, 1, function(x) {
            paste0(param, "[", paste(rev(x), collapse = ","), "]")
        }))
      } else
      {
        a <- list(y)
        names(a) <- param
      }

      ## loop over all parameters (if dimensionsless, just the parameter,
      ## otherwise all possible dimensions) and remove all dimension columns
      a <- lapply(names(a), function(x) {
        for (col in colnames(unique_dims)) {
          a[[x]][[col]] <- NULL
        }
        data.table::setnames(a[[x]], "value", x)
      })
      l <- c(l, a)
    }

    ## create a wide table of all the parameters, for calculating
    ## covariances
    wide <- l[[1]]
    if (length(l) > 1) {
      for (i in seq(2, length(l))) {
        wide <- merge(wide, l[[i]], by=intersect(colnames(wide), colnames(l[[i]])))
      }
    }
    wide[["np"]] <- NULL

    ## calculate the covariance matrix
    c <- stats::cov(wide)
    if (start) {
      c[, ] <- 0
    }

    ## calculate the vector of variances, and scaling of the mean
    sd_vec <- diag(c) - c[1, ]**2 / c[1, 1]
    mean_scale <- c[1, ] / c[1, 1]

    ## in case machine precision has made something < 0, set it to 0
    sd_vec[!(is.finite(sd_vec) & sd_vec > 0)] <- 0
    mean_scale[!(is.finite(mean_scale))] <- 0

    ## take square root of variances to get standard deviation
    sd_vec <- sqrt(sd_vec)
  } else {
    sd_vec <- vapply(all_params, function(p) {
      ifelse(length(res[[p]]) == 1, 0, sd(res[[p]]$value))
    }, 0)
  }

  if (missing(scale)) {
    scale_string <- ""
  } else {
    scale_string <- paste0(scale, " * ")
  }

  ## get prior density definition for each parameter
  param_bounds <-
    vapply(all_params, function(param) {
      grep(paste0("^[[:space:]]*", param, "[[[:space:]][^~]*~"),
           unlist(param_block), value = TRUE)
    }, "")
  ## select parameter that are variable (has a ~ in its prior line)
  variable_bounds <- param_bounds[vapply(param_bounds, function(x) {
    length(x) > 0
  }, TRUE)]

  proposal_lines <- list(parameter=c(), initial=c())

  first <- TRUE ## we're at the first parameter
  for (dim_param in names(sd_vec)) { ## loop over all parameters
    ## create dimensionless parameter
    param <- gsub("\\[[^]]*\\]", "", dim_param)
    ## check if it's variable
    if (param %in% names(variable_bounds[param]))
    {
      block <- names(params)[vapply(params, function(x) {param %in% x}, TRUE)]
      ## extract name of the parameter plus dimensions
      param_string <-
        sub(paste0("^[:space:]*(", param, "[^[:space:]~]*)[[:space:]~].*$"), "\\1", variable_bounds[param])
      ## extract bounded distribution split from parameters
      param_bounds_string <-
        sub("^.*(uniform|truncated_gaussian|truncated_normal|gamma|beta)\\((.*)\\)[[:space:]]*$",
            "\\1|\\2", variable_bounds[param])

      ## split distribution from arguments
      args <- strsplit(param_bounds_string, split = "\\|")
      ## extract distribution
      dist <- args[[1]][1]
      ## extract arguments to distribution
      bounds_string <- args[[1]][2]

      if (first) {
        ## first parameter; this is treated slightly different from the
        ## others in building up the multivariate normal distribution from
        ## interdependent univariate normal distributions
        mean <- ifelse(correlations, dim_param, param_string)
        if (correlations) {
          old_name <- "_old_mean_"
          proposal_lines[[block]] <- paste("inline", old_name, "=", dim_param)
          sd <- sqrt(c[dim_param, dim_param])
        } else {
          sd <- sd_vec[[dim_param]]
        }
      } else {
        if (correlations) {
          mean <- paste0(dim_param, " + (", mean_scale[dim_param], ") * __diff_")
        } else {
          mean <- param_string
        }
        sd <- sd_vec[[dim_param]]
      }

      ## impose bounds on gamma and beta distributions
      if (dist == "beta") {
        bounds_string <- "lower = 0, upper = 1"
      } else if (dist %in% c("gamma", "poisson", "negbin", "binomial")) {
        bounds_string <- "lower = 0"
      }

      if (is.na(bounds_string) || bounds_string == variable_bounds[param]) {
        ## no bounds, just use a gaussian
        if (sd == 0) {
          sd <- 1
        }
        proposal_lines[[block]] <-
          c(proposal_lines[[block]],
            paste0(ifelse(correlations, dim_param, param_string), " ~ gaussian(",
                   "mean = ", mean,
                   ", std = ", scale_string, sd, ")"))
      } else {
        ## there are (potentially) bounds, use a truncated normal
        bounds <- c(lower = NA, upper = NA)

        ## extract upper and lower bounds
        split_bounds <- strsplit(bounds_string, split = ",")[[1]]
        for (bound in c("lower", "upper")) {
          named <- grep(paste0(bound, "[[:space:]]*="), split_bounds)
          if (length(named) > 0) {
            bounds[bound] <- split_bounds[named]
            split_bounds <- split_bounds[-named]
          }
        }

        ## remove any arguments that don't pertain to bounds
        if (any(is.na(bounds))) {
          if (length(grep("^truncated", dist)) > 0) {
            named_other <- grep("(mean|std)", split_bounds)
            if (length(named_other) > 0) {
              split_bounds <- split_bounds[-named_other]
            }
            if (length(named_other) < 2) {
              split_bounds <- split_bounds[-seq_len(2 - length(named_other))]
            }
          }
        }

        ## get bounds
        for (split_bound in split_bounds) {
          bounds[which(is.na(bounds))][1] <- split_bound
        }

        ## get lower and upper bound
        bounds <- gsub("(lower|upper)[[:space:]]*=[[:space:]]*", "", bounds)
        bounds <- bounds[!is.na(bounds)]

        ## evaluate bounds (if they are given as expressions)
        eval_bounds <- tryCatch(
        {
          vapply(bounds, function(x) {as.character(eval(parse(text = x)))}, "")
        },
        error = function(cond)
        {
          ## preserve adapted dimensions
          if (param != dim_param) {
            orig_param_dims <- sub("^.*\\[(.*)\\]$", "\\1", param_string)
            orig_param_dims <- unlist(strsplit(orig_param_dims, ","))
            new_param_dims <- sub("^.*\\[(.*)\\]$", "\\1", dim_param)
            new_param_dims <- unlist(strsplit(new_param_dims, ","))
            names(new_param_dims) <- orig_param_dims
            for (bound in names(bounds))
            {
              orig_bound_dims <- sub("^.*\\[(.*)\\]$", "\\1", bounds[bound])
              bound_dims <- unlist(strsplit(orig_bound_dims, ","))
              matching_dims <- which(bound_dims %in% names(new_param_dims))
              if (length(matching_dims) > 0)
              {
                bound_dims[matching_dims] <-
                  new_param_dims[bound_dims[bound_dims %in% names(new_param_dims)]]
              }
              bounds[bound] <-
                sub(paste0("\\[", orig_bound_dims, "\\]"),
                    paste0("[",paste(bound_dims, collapse=","), "]"),
                    bounds[bound])
            }
            ret <- bounds
            return(ret)
          }
        })
        bounds <- eval_bounds

        if (sd == 0) {
          ## no variation
          if (sum(!is.na(as.numeric(bounds))) == 2) {
            ## range / 10
            sd <- diff(as.numeric(bounds)) / 10
          } else {
            sd <- 1
          }
        } else {
          if (sum(!is.na(as.numeric(bounds))) == 2) {
            sd <- min(sd, diff(as.numeric(bounds)) / 10)
          }
        }

        ## construct proposal
        proposal_lines[[block]] <-
          c(proposal_lines[[block]],
            paste0(ifelse(correlations, dim_param, param_string),
                   " ~ truncated_gaussian(",
                   "mean = ", mean,
                   ", std = ", scale_string, sd,
                   ifelse(length(bounds) > 0,
                          paste0(", ", paste(names(bounds), "=", bounds,
                                             sep = " ", collapse = ", "),
                                 ")"),
                          ")")))
      }

      if (first) {
        first <- FALSE
        if (correlations) {
          proposal_lines[[block]] <-
            c(proposal_lines[[block]],
              paste("__diff_ <- ", dim_param, "-", "_old_mean_"))
        }
      }
    }
  }

  for (block in c("parameter", "initial"))
  {
    model <- add_block(model, name = paste0("proposal_", block),
                       lines = proposal_lines[[block]])
  }

  return(model)
}
