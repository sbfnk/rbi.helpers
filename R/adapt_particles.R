#' @rdname adapt_particles
#' @name adapt_particles
#' @title Adapt the number of particles
#' @description This function takes the provided \code{\link{libbi}} and
#'   runs MCMC at a single point (i.e., repeatedly proposing the same paramters),
#'   adapting the number of particles distribution until the variance of the log-likelihood
#'   crosses the value given as \code{target.variance} (1 by default).
#' @param wrapper \code{\link{libbi}} (which has been run) to study
#' @param min minimum number of particles
#' @param max maximum number of particles
#' @param nsamples number of samples to generate each iteration; if not give, will use what has been set in \code{wrapper}
#' @param target.variance target log-likelihood variance; once this is crossed, the current number of particles will be used
#' @param ... parameters for libbi$run
#' @return a \code{\link{libbi}} with the desired proposal distribution
#' @importFrom coda mcmc rejectionRate effectiveSize
#' @importFrom rbi bi_model bi_dim_len bi_read
#' @importFrom stats var
#' @examples
#' example_obs <- bi_read(system.file(package="rbi", "example_output.nc"))
#' example_model <- bi_model(system.file(package="rbi", "PZ.bi"))
#' example_bi <- libbi(model = example_model, obs = example_obs)
#' obs_states <- example_model$get_vars("obs")
#' max_time <- max(sapply(example_obs[obs_states], function(x) { max(x[["time"]])}))
#' \dontrun{adapted <- adapt_particles(example_bi, nsamples = 128, end_time = max_time)}
#' @export
adapt_particles <- function(wrapper, min = 1, max = 1024, nsamples, target.variance = 1, quite=FALSE, ...) {

  if (!wrapper$run_flag) {
    message(date(), " Initial trial run")
    wrapper$run(client = "sample", ...)
  }

  if (missing(min) && !is.null(wrapper$options[["nparticles"]]))
  {
    min <- wrapper$options[["nparticles"]]
  }
  if (max <= min) {
    stop("'max' must be less or equal to 'min'")
  }

  test <- 2**(seq(floor(log(min, 2)), ceiling(log(max, 2))))

  model <- rbi::bi_model(lines = wrapper$model$get_lines())
  model$remove_block("proposal_parameter")
  model$remove_block("proposal_initial")
  model$remove_block("parameter")

  adapt_wrapper <- wrapper$clone(model = model)
  init_wrapper <- wrapper

  options <- list()
  ## use last parameter value from output file
  options[["init-np"]] <- rbi::bi_dim_len(wrapper$result$output_file_name, "np") - 1

  if (missing(nsamples)) {
    if ("nsamples" %in% names(adapt_wrapper$options)) {
      nsamples <- adapt_wrapper$options[["nsamples"]]
    } else {
      stop("if 'nsamples' is not an option in the 'libbi' object, it must be provided")
    }
  } else {
    options[["nsamples"]] <- nsamples
  }

  accRate <- c()
  var_loglik <- c()
  found_good <- FALSE
  id <- 0
  while (!found_good && id < length(test)) {
    id <- id + 1
    options[["nparticles"]] <- test[id]
    adapt_wrapper <-
      adapt_wrapper$clone(model = model, run = TRUE, options = options,
                          init = init_wrapper, ...)
    options[["init-np"]] <- nsamples - 1

    var_loglik <- c(var_loglik, stats::var(rbi::bi_read(adapt_wrapper, "loglikelihood")$value))

    message(date(), " ", test[id], " particles, loglikelihod variance: ", var_loglik[id])

    if (var_loglik[id] > 0) {
      init_wrapper <- adapt_wrapper
      if (var_loglik[id] < target.variance) {
      ## choose smallest var-loglikelihood < target.variance
        found_good <- TRUE
        if (id > 1) id <- id - 1
      }
    }
  }

  adapt_wrapper$options[["nparticles"]] <- test[id]
  adapt_wrapper$model <- wrapper$model

  message(date(), " Choosing ", test[id], " particles.")

  return(adapt_wrapper)
}
