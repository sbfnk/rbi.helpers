#' @rdname update_proposal
#' @name update_proposal
#' @title Construct a proposal from run results
#' @description
#' This function takes the provided \code{\link{bi_model}} and adds a generic
#' Gaussian proposal distribution.
#' @param model a \code{\link{bi_model}} object
#' @param correlations whether to take into account correlations
#' @param truncate truncate the multivariate normal proposals according to the used priors, e.g. truncating a parameter with beta prior at 0 and 1
#' @param blocks blocks to use (out of "parameter" and "initial")
#' @importFrom rbi get_block add_block insert_lines
#' @return the updated bi model
#' @keywords internal
update_proposal <- function(model, correlations = FALSE, truncate = TRUE, blocks=c("parameter", "initial")) {

  if (!("bi_model" %in% class(model))) stop("'x' must be a 'bi_model'")

  ## get constant expressions
  const <- get_const(model)
  for (const_name in names(const)) {
    assign(const_name, const[[const_name]])
  }

  param_block <- list()
  for (block in blocks)
  {
    ## get parameters
    param_block[[block]] <- get_block(model, block)
  }

  ## get dimensionless parameters
  params <- list(parameter=var_names(model, type="param"),
                 initial=var_names(model, type="state"))
  params <- params[blocks]

  ## get dimension parameters
  dims <- get_dims(model)
  block_to_state <- c(parameter="param", initial="state")
  vars <- var_names(model, type=block_to_state[blocks], dim=TRUE)
  dim_vars <- unlist(lapply(vars, function(x)
  {
    if (grepl("\\[", x)) {
      var_dims <-
        gsub("[[:space:]]", "", strsplit(sub("^.*\\[(.+)\\]$", "\\1", x), ",")[[1]])
      dim_list <- lapply(var_dims, function(x) seq_len(dims[[x]])-1)
      names(dim_list) <- var_dims

      df <- do.call(expand.grid, dim_list)
      all_dims <- apply(df, 1, function(x) paste(x, sep=","))
      var_name <- gsub("[[:space:]]", "", sub("\\[.+$", "", x))
      paste0(var_name, "[", all_dims, "]")
    } else {
      x
    }
  }))

  ## get prior density definition for each parameter
  param_bounds <-
    vapply(dim_vars, function(param) {
      match <- grep(paste0("^[[:space:]]*", gsub("([][])", "\\\\\\1", param),
                           "[[:space:]][^~]*~"),
                    unlist(param_block), value = TRUE)
      if (length(match) == 0) {
        param_regex <- gsub("\\[[0-9, ]*\\]", "[[A-z0-9_,: ].*]", param)
        match <- grep(paste0("^[[:space:]]*", param_regex, "[[:space:]][^~]*~"),
                      unlist(param_block), value = TRUE)
        if (length(match) > 1) {
          id <- as.integer(sub("^.*\\[([0-9, ]*)\\].*$", "\\1", param))
          ranges <- lapply(match, function(x) {
            eval(parse(text=sub("^.*\\[([ 0-9:,]*)\\].*~.*$", "\\1", x)))
          })
          match <- match[which(vapply(ranges, function(x) id %in% x, TRUE))]
        }
      }
      if (length(match) == 0) return("")
      return(match)
    }, "")
  param_bounds <- param_bounds[nchar(param_bounds) > 0]
  ## select parameter that are variable (has a ~ in its prior line)
  variable_bounds <- param_bounds[vapply(param_bounds, function(x) {
    length(x) > 0
  }, TRUE)]

  variable_names <- names(variable_bounds)
  dimless_variable_names <- sub("[[:space:]]*\\[.*$", "", variable_names)

  block_vars <- lapply(names(params), function(x)
  {
    variable_names[dimless_variable_names %in% params[[x]]]
  })
  dimless_block_vars <- lapply(names(params), function(x)
  {
    dimless_variable_names[dimless_variable_names %in% params[[x]]]
  })
  names(block_vars) <- names(params)
  names(dimless_block_vars) <- names(params)

  for (block in names(block_vars))
  {
    proposal_lines <- c()
    dim_lines <- grep("^[[:space:]]*dim[[:space:]]", model)
    if (length(dim_lines)==0) dim_lines <- 1

    for (param_id in seq_along(block_vars[[block]]))
    {
      param <- block_vars[[block]][param_id]
      dimless_param <- dimless_block_vars[[block]][param_id]
      dim_param <- var_names(model, vars=dimless_param, dim=TRUE)

      ## extract bounded distribution split from parameters
      param_bounds_string <-
        sub("^.*(uniform|truncated_gaussian|truncated_normal|gamma|beta|exponential|binomial|negbin|betabin)\\((.*)\\)[[:space:]]*$",
            "\\1|\\2", variable_bounds[param])

      ## split distribution from arguments
      args <- strsplit(param_bounds_string, split = "\\|")
      ## extract distribution
      dist <- args[[1]][1]
      ## extract arguments to distribution
      bounds_string <- args[[1]][2]

      mean <- param

      if (correlations) {
        for (j in (seq_len(param_id - 1)))
        {
          mean <-
            paste0(mean, " + __proposal_", block, "_cov[",
                   param_id-1, ",", j-1, "] * (",
                   block_vars[[block]][j], " - ",
                   paste0("__current_", block_vars[[block]][j]), ")")
        }
        std <- paste0("__proposal_", block, "_cov[", param_id-1, ",", param_id-1, "]")
      } else {
        mean <- param
        std <- paste0("__std_", param)
      }

      ## impose bounds on gamma and beta distributions
      if (dist == "beta") {
        bounds_string <- "lower = 0, upper = 1"
      } else if (dist %in% c("gamma", "poisson", "negbin", "binomial", "betabin", "exponential")) {
        bounds_string <- "lower = 0"
      }

      if (!truncate || is.na(bounds_string) ||
            bounds_string == variable_bounds[param]) {
        ## no bounds, just use a gaussian
        proposal_lines <-
          c(proposal_lines,
            paste0(param, " ~ gaussian(", "mean = ", mean, ", std = ", std, ")"))
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
          if (dimless_param != param) {
            orig_param_dims <- sub("^.*\\[(.*)\\]$", "\\1", dim_param)
            orig_param_dims <- unlist(strsplit(orig_param_dims, ","))
            new_param_dims <- sub("^.*\\[(.*)\\]$", "\\1", param)
            new_param_dims <- unlist(strsplit(new_param_dims, ","))
            names(new_param_dims) <- orig_param_dims
            for (bound in names(bounds))
            {
              orig_bound_dims <- sub("^.*\\[(.*)\\]$", "\\1", bounds[bound])
              bound_dims <- unlist(strsplit(orig_bound_dims, ","))
              matching_dims <- which(bound_dims %in% orig_param_dims)
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

        proposal_lines <-
          c(proposal_lines,
            paste0(param, " ~ truncated_gaussian(", "mean = ", mean, ", std = ", std,
                   ifelse(length(bounds) > 0,
                          paste0(", ", paste(names(bounds), "=", bounds,
                                             sep = " ", collapse = ", "),
                                 ")"),
                          ")")))
      }
    }

    var_type <- block_to_state[block]
    vars <- var_names(model, type=var_type, dim=FALSE)
    dim_vars <- var_names(model, type=var_type, dim=TRUE)
    propose_parameters <- rev(dim_vars[vars %in% unique(dimless_variable_names)])

    if (length(propose_parameters) > 0) {
      if (correlations) {
        proposal_lines <-
          c(paste0("__current_", propose_parameters, " <- ", propose_parameters),
            proposal_lines)
        new_param_names <-
          setdiff(paste0("__current_", propose_parameters), vars)
        if (length(new_param_names) > 0) {
          var_lines <- paste(var_type, new_param_names, "(has_output=0)")
          model <- insert_lines(model, var_lines, after=max(dim_lines))
        }
        new_dim <- paste0("__dim_", block, "_cov")
        cov_lines <-
          c(paste0("dim ", new_dim, "(", length(block_vars[[block]]), ")"),
            paste0("input __proposal_", block, "_cov[", new_dim, ",", new_dim, "]"))
        model <- insert_lines(model, cov_lines, after=max(dim_lines))
      } else {
        new_param_names <-
          paste0("__std_", propose_parameters)
        model <- insert_lines(model, paste("input", new_param_names),
                              after=max(dim_lines))
      }
    }

    model <- add_block(model, name = paste0("proposal_", block),
                       lines = proposal_lines)
  }

  return(model)
}
