#' @rdname adapt_particles
#' @name adapt_particles
#' @title Adapt the number of particles
#' @description This function takes the provided [rbi::libbi()] and runs
#'   MCMC at a single point (i.e., repeatedly proposing the same parameters),
#'   adapting the number of particles distribution until the variance of the
#'   log-likelihood crosses the value given as \code{target.variance} (1 by
#'   default).
#' @param x a [rbi::libbi()] object
#' @param min minimum number of particles
#' @param max maximum number of particles
#' @param target_variance target log-likelihood variance; once this is crossed,
#'   the current number of particles will be used
#' @param quiet if set to TRUE, will not provide running output of particle
#'   numbers tested
#' @param target.variance deprecated; use \code{target_variance} instead
#' @param ... parameters for libbi$run
#' @return a [rbi::libbi()] with the desired proposal distribution
#' @importFrom rbi bi_model bi_read sample remove_lines attach_data
#'   installed_libbi_version bi_contents var_names
#' @importFrom stats var
#' @importFrom utils compareVersion
#' @examples
#' example_obs <- rbi::bi_read(system.file(package="rbi", "example_dataset.nc"))
#' example_model <- rbi::bi_model(system.file(package="rbi", "PZ.bi"))
#' example_bi <- rbi::libbi(model = example_model, obs = example_obs)
#' obs_states <- rbi::var_names(example_model, type = "obs")
#' max_time <- max(vapply(example_obs[obs_states], function(x) {
#'   max(x[["time"]])
#' }, 0))
#' \dontrun{
#'   adapted <- adapt_particles(example_bi, nsamples = 128, end_time = max_time)
#' }
#' @export
adapt_particles <- function(x, min = 1, max = 1024, target_variance = 1,
                            quiet = FALSE, target.variance, ...) { # nolint
  if (!missing(target.variance)) {
    warning("target.variance is deprecated. Use 'target_variance' instead.")
    if (missing(target_variance)) {
      target_variance <- target.variance
    } else {
      stop("Can't give 'target.variance' and 'target_variance'.")
    }
  }

  if (missing(min) && !is.null(x$options[["nparticles"]])) {
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

  libbi_version <- installed_libbi_version(x$path_to_libbi)
  old_libbi <- ## check if version is > 1.4.2
    (compareVersion(libbi_version, "1.4.2") <= 0)

  thin <- x$thin ## no thinning when adapting particles

  test <- union(min * 2**(seq(0, floor(log2(max / min)))), max)

  adapted <- x

  if (old_libbi) { ## old libbi version; we must copy the prior to the proposal
    model <- adapted$model
    model <- rbi::remove_lines(model, "proposal_parameter")
    model <- rbi::remove_lines(model, "parameter")
    if ("with-transform-initial-to-param" %in% names(x$options)) {
      model <- rbi::remove_lines(model, "proposal_initial")
      model <- rbi::remove_lines(model, "initial")
    }
    adapted$model <- model
  } else {
    adapted$model <- update_proposal(adapted$model)
    zero_cov_vars <- get_mvn_params(adapted, fix = 0)
    if ("input-file" %in% names(adapted$options)) {
      existing_cov_vars <-
        intersect(bi_contents(adapted, file = "input"), names(zero_cov_vars))
      existing_cov_input <- bi_read(
        adapted, file = "input", vars = existing_cov_vars
      )
    } else {
      existing_cov_input <- list()
    }
    adapted <- attach_data(
      adapted, file = "input", zero_cov_vars, append = TRUE, overwrite = TRUE
    )
  }

  adapted$thin <- 1

  found_good <- FALSE
  id <- 0
  while (!found_good && (id + 1) < length(test)) {
    id <- id + 1
    adapted <- rbi::sample(adapted, nparticles = test[id], chain = TRUE, ...)

    loglik <- rbi::bi_read(adapted, "loglikelihood")$loglikelihood$value
    ## remove any infinite log-likelihoods
    loglik <- loglik[is.finite(loglik)]
    var_loglik <- stats::var(loglik)
    if (is.na(var_loglik)) var_loglik <- 0

    if (!quiet) {
      message(
        date(), " ", test[id], " particles, loglikelihod variance: ",
        var_loglik
      )
    }

    if (var_loglik < target_variance) {
      ## choose smallest var-loglikelihood < target_variance
      found_good <- TRUE
    }
  }
  if (!found_good) id <- length(test)

  adapted$options[["nparticles"]] <- test[id]
  adapted$thin <- thin

  if (old_libbi) {
    adapted$model <- x$model
    adapted$run_flag <- FALSE
  } else if (length(existing_cov_input) > 0) {
    adapted <-
      attach_data(adapted,
        file = "input", existing_cov_input,
        append = TRUE, overwrite = TRUE
      )
  }

  if (!quiet) message(date(), " Choosing ", test[id], " particles.")

  return(adapted)
}
