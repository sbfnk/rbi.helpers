#' @rdname adapt_particles
#' @name adapt_particles
#' @title Adapt the number of particles 
#' @description This function takes the provided \code{\link{libbi}} and
#'   runs MCMC at a single point (i.e., repeatedly proposing the same paramters),
#'   adapting the number of particles distribution until the desired 
#'   acceptance rate is achieved. If a scale is given, it will be used
#'   to adapt the number of particles at each iteration.
#' @param wrapper \code{\link{libbi}} (which has been run) to study
#' @param add_options list of additional options
#' @param samples number of samples to generate each iteration
#' @param ... parameters for libbi$run
#' @importFrom coda mcmc rejectionRate effectiveSize
#' @return a \code{\link{libbi}} with the desired proposal distribution
#' @export
adapt_particles <- function(wrapper, test = 2**(0:10), add_options, samples, ...) {

  if (missing(add_options)) {
    add_options <- list()
  } else if (!is.list(add_options)) {
    stop("'add_options' must be given as list.")
  }

  if (!wrapper$run_flag) {
    stop("The model should be run first")
  }

  model <- bi_model(lines = wrapper$model$get_lines())
  model$remove_block("proposal_parameter")
  model$remove_block("proposal_initial")
  model$remove_block("parameter")

  adapt_wrapper <- wrapper$clone(model = model)
  init_wrapper <- wrapper

  ## use last parameter value from output file
  add_options[["init-np"]] <- bi_dim_len(wrapper$output_file_name, "np") - 1

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
  ess <- c()
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
    init_wrapper <- adapt_wrapper
    
    mcmc_obj <- mcmc(get_traces(adapt_wrapper, all = TRUE))
    accRate <- c(accRate, max(1 - rejectionRate(mcmc_obj)))
    ess <- c(ess, max(effectiveSize(mcmc_obj)))
    var_loglik <- c(var_loglik, var(bi_read(adapt_wrapper, "loglikelihood")$value))
    
    cat(paste0(test[id], " particles: acceptance rate ", accRate[id],
               ", ESS ", ess[id], ", loglikelihod variance: ", var_loglik[id], "\n"))

    ## choose smallest var-loglikelihood < 1
    if (var_loglik[id] < 1) {
      found_good <- TRUE
      if (id > 1) id <- id - 1
    }
  }

  adapt_wrapper$global_options[["nparticles"]] <- test[id]
  adapt_wrapper$model <- wrapper$model

  cat("Choosing ", test[id], " particles.\n")
      
  return(adapt_wrapper)
}
