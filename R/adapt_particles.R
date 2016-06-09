#' @rdname adapt_particles
#' @name adapt_particles
#' @title Adapt the number of particles 
#' @description This function takes the provided \code{\link{libbi}} and
#'   runs MCMC at a single point (i.e., repeatedly proposing the same paramters),
#'   adapting the number of particles distribution until the variance of the likelihood
#'   crosses one.
#' @param wrapper \code{\link{libbi}} (which has been run) to study
#' @param min minimum number of particles
#' @param max maximum number of particles
#' @param add_options list of additional options
#' @param samples number of samples to generate each iteration
#' @param ... parameters for libbi$run
#' @importFrom coda mcmc rejectionRate effectiveSize
#' @return a \code{\link{libbi}} with the desired proposal distribution
#' @export
adapt_particles <- function(wrapper, min = 1, max = 1024, add_options, samples, ...) {

  if (missing(add_options)) {
    add_options <- list()
  } else if (!is.list(add_options)) {
    stop("'add_options' must be given as list.")
  }

  if (!wrapper$run_flag) {
    stop("The model should be run first")
  }

  if (missing(min) && !is.null(wrapper$global_options[["nparticles"]]))
  {
    min <- wrapper$global_options[["nparticles"]]
  }
  if (max <= min) {
    stop("'max' must be less or equal to 'min'")
  }

  test <- 2**(seq(floor(log(min, 2)), ceiling(log(max, 2))))

  model <- bi_model(lines = wrapper$model$get_lines())
  model$remove_block("proposal_parameter")
  model$remove_block("proposal_initial")
  model$remove_block("parameter")

  adapt_wrapper <- wrapper$clone(model = model)
  init_wrapper <- wrapper

  ## use last parameter value from output file
  add_options[["init-np"]] <- bi_dim_len(wrapper$result$output_file_name, "np") - 1

  if (missing(samples)) {
    if ("nsamples" %in% names(adapt_wrapper$global_options)) {
      samples <- adapt_wrapper$global_options[["nsamples"]]
    } else {
      stop("if 'nsamples' is not a global option, must provide 'samples'")
    }
  } else {
    add_options[["nsamples"]] <- samples
  }

  accRate <- c()
  var_loglik <- c()
  found_good <- FALSE
  id <- 0
  while (!found_good && id < length(test)) {
    id <- id + 1
    add_options[["nparticles"]] <- test[id]
    adapt_wrapper <-
      adapt_wrapper$clone(model = model, run = TRUE, add_options = add_options,
                          init = init_wrapper, ...)
    add_options[["init-np"]] <- samples - 1

    var_loglik <- c(var_loglik, var(bi_read(adapt_wrapper, "loglikelihood")$value))

    cat(date(), paste0(test[id], ", loglikelihod variance: ", var_loglik[id], "\n"))

    if (var_loglik[id] > 0) {
      init_wrapper <- adapt_wrapper
      if (var_loglik[id] < 1) {
      ## choose smallest var-loglikelihood < 1
        found_good <- TRUE
        if (id > 1) id <- id - 1
      }
    }
  }

  adapt_wrapper$global_options[["nparticles"]] <- test[id]
  adapt_wrapper$model <- wrapper$model

  cat(date(), "Choosing ", test[id], " particles.\n")

  return(adapt_wrapper)
}
