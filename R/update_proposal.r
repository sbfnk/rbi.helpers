#' @rdname update_proposal
#' @name update_proposal
#' @title Construct a proposal from run results
#' @description
#' This function takes the provided \code{\link{bi_model}} and adds a generic
#' Gaussian proposal distribution.
#' @param model a \code{\link{bi_model}} object
#' @param truncate truncate the multivariate normal proposals according to the used priors, e.g. truncating a parameter with beta prior at 0 and 1
#' @param blocks blocks to use (out of "parameter" and "initial")
#' @importFrom rbi get_block add_block insert_lines get_const rewrite
#' @importFrom utils capture.output
#' @return the updated bi model
#' @keywords internal
update_proposal <- function(model, truncate = TRUE, blocks=c("parameter", "initial")) {

  if (!("bi_model" %in% class(model))) stop("'model' must be a 'bi_model'")

  ## get constant expressions
  const <- get_const(model)
  for (const_name in names(const)) {
    assign(const_name, const[[const_name]])
  }

  libbi_model <- bi_model(lines=capture.output(rewrite(model)))
  param_block <- list()
  ranges <- list()

  for (block in blocks)
  {
    ## get parameters
    param_block[[block]] <- grep("~", get_block(libbi_model, block),  value=TRUE)
    ranges[[block]] <- unlist(lapply(param_block[[block]], function(x)
    {
      ## strip to parameter definition
      dim_param <- sub("[[:space:]]*~.*$", "", x)
      ## strip to distribution
      dist_param <- sub("^.*~[[:space:]]*", "", x)
      ## strip dim identifiers
      dim_param <- sub("[^[,]+=", "", dim_param)
      ret <- list(param=sub("\\[.*$", "", dim_param), dist=dist_param)
      if (grepl("\\[", dim_param)) {
        vardims <- unlist(strsplit(sub("^.*\\[", "", sub("\\]$", "", dim_param)), ","))
        grid <- do.call(expand.grid, lapply(vardims, function(x) eval(parse(text=x))))
        ret <- apply(grid, 1, function(x)
        {
          c(ret, list(dim=paste0("[", paste(x, collapse=","), "]")))
        })
      } else {
        ret <- list(c(ret, list(dim="")))
      }
      ret
    }), recursive = FALSE)
  }

  for (block in names(ranges))
  {
    proposal_lines <- c()
    dim_lines <- grep("^[[:space:]]*dim[[:space:]]", model)
    if (length(dim_lines)==0) dim_lines <- 1
    block_vars <- c()

    for (param_id in seq_along(ranges[[block]]))
    {
      range <- ranges[[block]][[param_id]]
      param <- paste0(range$param, range$dim)
      block_vars <- c(block_vars, param)
      dim_param <- var_names(model, vars=range$param, dim=TRUE)

      ## extract bounded distribution split from parameters
      param_bounds_string <-
        sub(paste0("^(uniform|truncated_gaussian|truncated_normal|gamma|beta|",
                   "exponential|binomial|negbin|betabin)\\((.*)\\)$"), "\\1|\\2",
            range$dist)

      ## split distribution from arguments
      args <- strsplit(param_bounds_string, split = "\\|")
      ## extract distribution
      dist <- args[[1]][1]
      ## extract arguments to distribution
      arg_string <- args[[1]][2]

      mean <- param

      for (j in (seq_len(param_id - 1))) {
        mean <-
          paste0(mean, " + __proposal_", block, "_cov[",
                 param_id-1, ",", j-1, "] * (", block_vars[j], " - ",
                 paste0("__current_", block_vars[j]), ")")
      }
      std <- paste0("__proposal_", block, "_cov[", param_id-1, ",", param_id-1, "]")

      ## impose bounds on gamma and beta distributions
      if (dist == "beta") {
        arg_string <- "lower = 0, upper = 1"
      } else if (dist %in% c("gamma", "poisson", "negbin", "binomial", "betabin", "exponential")) {
        arg_string <- "lower = 0"
      }

      ## extract upper and lower bounds (ignoring commas inside brackets)
      ## see https://stackoverflow.com/questions/39733645/split-string-by-space-except-whats-inside-parentheses
      split_args <-
        unlist(strsplit(arg_string,
                        "(\\[(?:[^[\\]]++|(?1))*\\])(*SKIP)(*F)|, ", perl=TRUE))
      named_args <- strsplit(split_args, split = " = ")
      arg_names <- vapply(named_args, function(x) x[[1]], "")
      bound_arg_pos <- which(arg_names %in% c("lower", "upper"))

      bound_args <- named_args[bound_arg_pos]

      if (!truncate || length(bound_args) == 0) {
        ## no bounds, just use a gaussian
        proposal_lines <-
          c(proposal_lines,
            paste0(param, " ~ gaussian(", "mean = ", mean, ", std = ", std, ")"))
      } else {
        ## evaluate bounds (if they are given as expressions)
        eval_bounds <- tryCatch(
        {
          bounds <-
            vapply(bound_args, function(x) as.character(eval(parse(text = x[[2]]))), "")
          names(bounds) <- vapply(bound_args, function(x) x[[1]], "")
          bounds
        },
        error = function(cond)
        {
          ## preserve adapted dimensions
          if (range$param != param) {
            orig_param_dims <- sub("^.*\\[(.*)\\]$", "\\1", dim_param)
            orig_param_dims <- unlist(strsplit(orig_param_dims, ","))
            new_param_dims <- sub("^.*\\[(.*)\\]$", "\\1", param)
            new_param_dims <- unlist(strsplit(new_param_dims, ","))
            names(new_param_dims) <- orig_param_dims
            bounds <- c()
            for (bound_id in seq_along(bound_args))
            {
              bound <- bound_args[[bound_id]][1]
              value <- bound_args[[bound_id]][2]
              orig_bound_dims <- sub("^.*\\[(.*)\\]$", "\\1", value)
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
                    value)
            }
            return(bounds)
          }
        })

        proposal_lines <-
          c(proposal_lines,
            paste0(param, " ~ truncated_gaussian(", "mean = ", mean, ", std = ", std,
                   paste0(", ", paste(names(eval_bounds), "=", eval_bounds,
                                      sep = " ", collapse = ", ")),
                   ")")
            )
      }
    }

    dims <- get_dims(model)
    block_to_state <- c(parameter="param", initial="state")
    var_type <- block_to_state[block]
    vars <- var_names(model, type=var_type, dim=FALSE)
    aux_inputs <- var_names(model, type="input", dim=FALSE, aux=TRUE)
    dim_vars <- var_names(model, type=var_type, dim=TRUE)
    dim_aux_vars <- var_names(model, type=var_type, dim=TRUE, aux=TRUE)
    dimless_variable_names <- unique(vapply(ranges[[block]], function(x) x$param, ""))
    propose_parameters <- rev(dim_vars[vars %in% unique(dimless_variable_names)])

    if (length(propose_parameters) > 0) {
      proposal_lines <-
        c(paste0("__current_", propose_parameters, " <- ", propose_parameters),
          proposal_lines)

      new_param_names <-
        setdiff(paste0("__current_", propose_parameters), dim_aux_vars)
      if (length(new_param_names) > 0) {
        new_var_lines <- paste(var_type, new_param_names, "(has_output=0)")
        model <- insert_lines(model, new_var_lines, after=max(dim_lines))
      }
      new_input_names <- setdiff(paste0("__proposal_", block, "_cov"), aux_inputs)
      cov_dim_name <- paste0("__dim_", block, "_cov")
      if (length(new_input_names) > 0) {
        new_input_lines <- paste0("input ", new_input_names, "[", cov_dim_name, ",", cov_dim_name, "]")
        model <- insert_lines(model, new_input_lines, after=max(dim_lines))
      }
      new_dim_names <- setdiff(cov_dim_name, names(dims))
      if (length(new_dim_names) > 0) {
        new_dim_lines <- paste0("dim ", new_dim_names, "(", length(block_vars), ")")
        model <- insert_lines(model, new_dim_lines, after=max(dim_lines))
      }
    }

    model <- add_block(model, name = paste0("proposal_", block),
                       lines = proposal_lines)
  }

  return(model)
}
