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
#' @param truncate truncate the multivariate normal proposals
#' @importFrom data.table setnames
#' @importFrom stats cov
#' @importFrom rbi get_block add_block insert_lines
#' @return the updated bi model
#' @keywords internal
output_to_proposal <- function(wrapper, scale, correlations = FALSE, truncate = TRUE) {

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

  param_block <- list()
  params <- list(parameter=var_names(model, "param"), initial=var_names(model, "state"))
  for (block in c("parameter", "initial"))
  {
    ## get parameters
    param_block[[block]] <- get_block(model, block)
  }

  ## read parameters
  init_params <- ("with-transform-initial-to-param" %in% names(wrapper$options))
  res_types <- c("param")
  if (init_params) res_types <- c(res_types, "state")
  res <- bi_read(wrapper, type=res_types, init.to.param = init_params)

  l <- list()
  for (param in names(res))
  {
    if (prod(dim(res[[param]])) > 0)
    {
      y <- data.table(res[[param]])
      ## extract columns that are dimensions
      dim_colnames <- setdiff(colnames(y), c("np", "value"))
      unique_dims <- unique(y[, dim_colnames, with=FALSE])
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
        if (length(a) > 0)
        {
          names(a) <- unname(apply(unique_dims, 1, function(x) {
            paste0(param, "[", paste(x, collapse = ","), "]")
          }))
        }
      } else
      {
        a <- list(y)
        names(a) <- param
      }
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

  ## get prior density definition for each parameter
  param_bounds <-
    vapply(colnames(wide), function(param) {
      match <- grep(paste0("^[[:space:]]*", gsub("([][])", "\\\\\\1", param), "[[:space:]][^~]*~"),
                    unlist(param_block), value = TRUE)
      if (length(match) == 0) {
        param_regex <- gsub("\\[[0-9, ]*\\]", "[[A-z_, ].*]", param)
        match <- grep(paste0("^[[:space:]]*", param_regex, "[[:space:]][^~]*~"),
                      unlist(param_block), value = TRUE)
      }
      if (length(match) == 0) return("")
      return(match)
    }, "")
  param_bounds <- param_bounds[nchar(param_bounds) > 0]
  ## select parameter that are variable (has a ~ in its prior line)
  variable_bounds <- param_bounds[vapply(param_bounds, function(x) {
    length(x) > 0
  }, TRUE)]

  if (length(variable_bounds) == 0) stop("No varying parameters. Cannot create proposal.")

  block_vars <- lapply(names(params), function(x)
  {
    variable_names <- names(variable_bounds)
    bare_variable_names <- sub(" *\\[.*$", "", variable_names)
    variable_names[bare_variable_names %in% params[[x]]]
  })
  names(block_vars) <- names(params)

  for (block in names(block_vars))
  {
    wide_block <- wide[, block_vars[[block]], with=FALSE]

    if (ncol(wide_block) > 0)
    {

      if (ncol(wide_block) > 1) {
        ## calculate the covariance matrix
        c <- stats::cov(wide_block)

        cholStatus <- try(A <- t(chol(c)), silent = TRUE)
        if (class(cholStatus) == "try-error")
        {
          warning("Cholesky decomposition failed; ",
                  "will try to adapt via independent univariate sampling first.")
          A <- NULL
        } else
        {
          A <- 2.38 * A / sqrt(ncol(A))
          Adiag <- diag(A)
          Afactors <- apply(A, 2, function(x) {x/Adiag})
        }
      } else {
        A <- NULL
      }

      sd_vec <- apply(wide_block, 2, sd)
      if (all(sd_vec == 0))
      {
        param_names <- names(sd_vec)
      } else {
        param_names <- names(sd_vec[sd_vec > 0])
      }
      dimless_param_names <- gsub("\\[[^]]*\\]", "", param_names)

      if (missing(scale)) {
        scale_string <- ""
      } else {
        scale_string <- paste0(scale, " * ")
      }

      proposal_lines <- c()

      for (param_id in seq_along(param_names)) { ## loop over all parameters
        dim_param <- param_names[param_id]
        param <- dimless_param_names[param_id]

        ## extract bounded distribution split from parameters
        param_bounds_string <-
          sub("^.*(uniform|truncated_gaussian|truncated_normal|gamma|beta)\\((.*)\\)[[:space:]]*$",
              "\\1|\\2", variable_bounds[dim_param])

        ## split distribution from arguments
        args <- strsplit(param_bounds_string, split = "\\|")
        ## extract distribution
        dist <- args[[1]][1]
        ## extract arguments to distribution
        bounds_string <- args[[1]][2]

        mean <- dim_param

        if (correlations && !is.null(A)) {
          for (j in seq_len(param_id - 1))
          {
            mean <- paste0(mean, " + (", Afactors[param_id, j], ") * (",
                           param_names[j], " - ", paste0("__", param_names[j]),
                           ")")
          }
          sd <- Adiag[[param_id]]
        } else {
          mean <- dim_param
          sd <- sd_vec[[dim_param]]
        }

        ## impose bounds on gamma and beta distributions
        if (dist == "beta") {
          bounds_string <- "lower = 0, upper = 1"
        } else if (dist %in% c("gamma", "poisson", "negbin", "binomial")) {
          bounds_string <- "lower = 0"
        }

        if (!truncate || is.na(bounds_string) || bounds_string == variable_bounds[dim_param]) {
          if (sd == 0) {
            sd <- 1
          }

          ## no bounds, just use a gaussian
          proposal_lines <-
            c(proposal_lines,
              paste0(dim_param, " ~ gaussian(", "mean = ", mean, ", std = ", scale_string, sd, ")"))
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
          }

          proposal_lines <-
            c(proposal_lines,
              paste0(dim_param, " ~ truncated_gaussian(", "mean = ", mean,
                     ", std = ", scale_string, sd,
                     ifelse(length(bounds) > 0,
                            paste0(", ", paste(names(bounds), "=", bounds,
                                               sep = " ", collapse = ", "),
                                   ")"),
                            ")")))
        }
      }

      ## HERE: use "param"/"states", not just "param"
      if (correlations && !is.null(A))
      {
        var_type <- ifelse(block == "parameter", "param", "state")
        vars <- var_names(model, type=var_type, dim=FALSE)
        dim_vars <- var_names(model, type=var_type, dim=TRUE)
        for (loop_dim_param in rev(dim_vars[vars %in% unique(dimless_param_names)]))
        {
          proposal_lines <- c(paste0("__", loop_dim_param, " <- ", loop_dim_param), proposal_lines)
        }
        for (loop_dim_param in rev(dim_vars[vars %in% unique(dimless_param_names)]))
        {
          proposal_lines <-
            c(paste(var_type, paste0("__", loop_dim_param), "(has_output=0)"), proposal_lines)
        }
      }

      model <- add_block(model, name = paste0("proposal_", block),
                         lines = proposal_lines)
    }
  }

  return(model)
}
