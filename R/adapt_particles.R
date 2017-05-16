#' @rdname adapt_particles
#' @name adapt_particles
#' @title Adapt the number of particles
#' @description This function takes the provided \code{\link{libbi}} and
#'   runs MCMC at a single point (i.e., repeatedly proposing the same paramters),
#'   adapting the number of particles distribution until the variance of the log-likelihood
#'   crosses the value given as \code{target.variance} (1 by default).
#' @param x a \code{\link{libbi}} object
#' @param min minimum number of particles
#' @param max maximum number of particles
#' @param target.variance target log-likelihood variance; once this is crossed, the current number of particles will be used
#' @param quiet if set to TRUE, will not provide running output of particle numbers tested
#' @param ... parameters for libbi$run
#' @return a \code{\link{libbi}} with the desired proposal distribution
#' @importFrom coda mcmc rejectionRate effectiveSize
#' @importFrom rbi bi_model bi_dim_len bi_read sample remove_lines
#' @importFrom stats var
#' @examples
#' example_obs <- bi_read(system.file(package="rbi", "example_dataset.nc"))
#' example_model <- bi_model(system.file(package="rbi", "PZ.bi"))
#' example_bi <- libbi(model = example_model, obs = example_obs)
#' obs_states <- var_names(example_model, "obs")
#' max_time <- max(vapply(example_obs[obs_states], function(x) { max(x[["time"]])}, 0))
#' \dontrun{adapted <- adapt_particles(example_bi, nsamples = 128, end_time = max_time)}
#' @export
adapt_particles <- function(x, min = 1, max = 1024, target.variance = 1, quiet=FALSE, ...) {

  if (missing(min) && !is.null(x$options[["nparticles"]]))
  {
    min <- x$options[["nparticles"]]
  }
  if (max < min) {
    stop("'max' must be greater or equal to 'min'")
  }

  if (!quiet) message(date(), " Adapting the number of particles")

  if (!x$run_flag) {
    if (!quiet) message(date(), " Initial trial run")
    x <- rbi::sample(x, ...)
  }

  thin <- x$thin ## no thinning when adapting particles

  test <- 2**(seq(floor(log(min, 2)), ceiling(log(max, 2))))

  adapted <- x
  model <- x$model
  model <- rbi::remove_lines(model, "proposal_parameter")
  model <- rbi::remove_lines(model, "parameter")
  if ("with-transform-initial-to-param" %in% names(x$options)) {
    model <- rbi::remove_lines(model, "proposal_initial")
    model <- rbi::remove_lines(model, "initial")
  }

  adapted$model <- model
  adapted$thin <- 1

  accRate <- c()
  var_loglik <- c()

  found_good <- FALSE
  id <- 0
  while (!found_good && id < length(test)) {
    id <- id + 1
    adapted <- rbi::sample(adapted, nparticles=test[id], chain=TRUE, ...)

    var_loglik <- c(var_loglik, stats::var(rbi::bi_read(adapted, "loglikelihood")$loglikelihood$value))

    if (!quiet) message(date(), " ", test[id], " particles, loglikelihod variance: ", var_loglik[id])

    if (var_loglik[id] < target.variance) {
      ## choose smallest var-loglikelihood < target.variance
      found_good <- TRUE
    }
  }

  adapted$options[["nparticles"]] <- test[id]
  adapted$model <- x$model
  adapted$thin <- thin

  if (!quiet) message(date(), " Choosing ", test[id], " particles.")

  return(adapted)
}
