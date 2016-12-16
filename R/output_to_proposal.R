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
#' @return the updated bi model
output_to_proposal <- function(wrapper, scale, correlations = FALSE, start = FALSE) {

  if (!wrapper$run_flag) {
    stop("The model should be run first")
  }

  model <- wrapper$model$clone()
  ## get constant expressions
  const_lines <- grep("^[[:space:]]*const", model$get_lines(), value = TRUE)
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
      warning("Cannot convert const expression for ", assignemnt[1],
              "into R expression")
      warning("Original message: ", cond)
    })
  }
  ## get parameters
  param_block <- model$get_block("parameter")
  ## only go over variable parameters
  params <- sub("[[:space:]]*~.*$", "", grep("~", param_block, value = TRUE))
  ## remove any dimensions
  params <- gsub("[[:space:]]*\\[[^]]*\\]", "", params)
  ## read parameters
  res <- bi_read(wrapper$output_file_name, vars = params)

  if (correlations) {
    l <- list()
    for (param in names(res)) {
      y <- copy(res[[param]])
      unique_dims <- unique(y[setdiff(colnames(y), c("np", "value"))])
      if (sum(dim(unique_dims)) > 0)
      {
        a <- apply(unique_dims, 1, function(x) {
          merge(t(x), y)
        })
        if (length(a))
          names(a) <- unname(apply(unique_dims, 1, function(x) {
            paste0(param, "[", paste(rev(x), collapse = ","), "]")
          }))
      } else
      {
        a <- list(y)
        names(a) <- param
      }

      a <- lapply(names(a), function(x) {
        for (col in colnames(unique_dims)) {
          a[[x]][[col]] <- NULL
        }
        data.table::setnames(a[[x]], "value", x)
      })
      l <- c(l, a)
    }

    if (length(l) > 1) {
      wide <- l[[1]]
      for (i in seq(2, length(l))) {
        wide <- merge(wide, l[[i]])
      }
    } else {
      wide <- l
    }
    wide[["np"]] <- NULL

    C <- stats::cov(wide)
    if (start) {
      C[, ] <- 0
    }

    sd_vec <- diag(C) - C[1, ]**2 / C[1, 1]
    mean_scale <- C[1, ] / C[1, 1]

    ## in case machine precision has made something < 0:
    sd_vec[!(is.finite(sd_vec) & sd_vec > 0)] <- 0
    mean_scale[!(is.finite(mean_scale))] <- 0
    sd_vec <- sqrt(sd_vec)
  } else {
    sd_vec <- sapply(params, function(p) {
      ifelse(length(res[[p]]) == 1, 0, sd(res[[p]]$value))
    })
  }

  if (missing(scale)) {
    scale_string <- ""
  } else {
    scale_string <- paste0(scale, " * ")
  }

  param_bounds <- sapply(params, function(param) {grep(paste0("^[[:space:]]*", param, "[[[:space:]][^~]*~"), param_block, value = TRUE)})
  variable_bounds <- param_bounds[sapply(param_bounds, function(x) {length(x) > 0})]

  proposal_lines <- c()

  first <- TRUE
  for (dim_param in names(sd_vec)) {
    param <- gsub("\\[[^]]*\\]", "", dim_param)
    if (param %in% names(variable_bounds[param]))
    {
      param_string <-
        sub(paste0("^[:space:]*(", param, "[^[:space:]~]*)[[:space:]~].*$"), "\\1", variable_bounds[param])
      param_bounds_string <-
        sub("^.*(uniform|truncated_gaussian|truncated_normal)\\((.*)\\)[[:space:]]*$",
            "\\1|\\2", variable_bounds[param])

      args <- strsplit(param_bounds_string, split = "\\|")

      dist <- args[[1]][1]
      bounds_string <- args[[1]][2]

      if (first) {
        mean <- ifelse(correlations, dim_param, param_string)
        if (correlations) {
          old_name <- "_old_mean_"
          proposal_lines <- paste("inline", old_name, "=", dim_param)
          sd <- sqrt(C[dim_param, dim_param])
        } else {
          sd <- sd_vec[[dim_param]]
        }
      } else {
        if (correlations) {
          mean <- paste0(dim_param, " + (", mean_scale[dim_param], ") * _old_mean_diff_")
        } else {
          mean <- param_string
        }
        sd <- sd_vec[[dim_param]]
      }

      if (is.na(bounds_string) || bounds_string == variable_bounds[param]) {
        if (sd == 0) {
          sd <- 1
        }
        proposal_lines <-
          c(proposal_lines,
            paste0(ifelse(correlations, dim_param, param_string), " ~ gaussian(",
                   "mean = ", mean,
                   ", std = ", scale_string, sd, ")"))
      } else {
        bounds <- c(lower = NA, upper = NA)

        split_bounds <- strsplit(bounds_string, split = ",")[[1]]
        for (bound in c("lower", "upper")) {
          named <- grep(paste0(bound, "[[:space:]]*="), split_bounds)
          if (length(named) > 0) {
            bounds[bound] <- split_bounds[named]
            split_bounds <- split_bounds[-named]
          }
        }

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

        for (split_bound in split_bounds) {
          bounds[which(is.na(bounds))][1] <- split_bound
        }

        bounds <- gsub("(lower|upper)[[:space:]]*=[[:space:]]*", "", bounds)
        bounds <- bounds[!is.na(bounds)]

        eval_bounds <- tryCatch(
        {
          sapply(bounds, function(x) { eval(parse(text = x))})
        },
        error = function(cond)
        {
          warning("Cannot convert bounds for ", param, "into R expression")
          warning("Original message: ", cond)
          ret <- bounds
          ret[] <- NA
          return(ret)
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
            sd <- min(sd, diff(as.numeric(bounds)) / 10
                      )
          }
        }

        proposal_lines <- c(proposal_lines,
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

      if (first && correlations) {
        proposal_lines <- c(proposal_lines, paste("inline", "_old_mean_diff_", "=", dim_param, "-", "_old_mean_"))
        first <- FALSE
      }
    }
  }

  model$add_block(name = "proposal_parameter",
                  lines = proposal_lines)

  return(model)
}
