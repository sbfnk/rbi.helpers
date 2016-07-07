#' @rdname adapt_mcmc
#' @name adapt_mcmc
#' @title Adapt the proposal distribution of MCMC using the covariance
#'   of samples
#' @description This function takes the provided \code{\link{libbi}} object and
#'   runs MCMC, adapting the proposal distribution until the desired
#'   acceptance rate is achieved. If a scale is given, it will be used
#'   to adapt the proposal at each iteration
#' @param wrapper \code{link{libbi}} (which has been run) to study
#' @param min minimum acceptance rate
#' @param max maximum acceptance rate
#' @param scale scale multiplier/divider for the proposal. If >1 this
#'   will be inverted.
#' @param add_options list of additional options
#' @param samples number of samples to generate each iteration
#' @param max_iter maximum of iterations (default: 10)
#' @param correlations take into account correlations
#' @param ... parameters for libbi$run
#' @return a \code{\link{libbi}} with the desired proposal distribution
#' @export
adapt_mcmc <- function(wrapper, min = 0, max = 1, scale = 2, add_options, samples, max_iter = 10, correlations = FALSE, ...) {

  if (missing(add_options)) {
    add_options <- list()
  } else if (!is.list(add_options)) {
    stop("'add_options' must be given as list.")
  }

  if (!wrapper$run_flag) {
    stop("The model should be run first")
  }

  ## scale should be > 1 (it's a divider if acceptance rate is too
  ## small, multiplier if the acceptance Rate is too big)
  if (scale < 1) scale <- 1 / scale

  init_file <- wrapper$result$output_file_name
  init_np <- bi_dim_len(init_file, "np") - 1

  if (missing(samples)) {
    if ("nsamples" %in% names(wrapper$global_options)) {
      samples <- wrapper$global_options[["nsamples"]]
    } else {
      stop("if 'nsamples' is not a global option, must provide 'samples'")
    }
  } else {
    add_options[["nsamples"]] <- samples
  }
  add_options[["init-file"]] <- init_file
  add_options[["init-np"]] <- init_np
  accRate <- acceptance_rate(wrapper)
  adapt_wrapper <- wrapper
  shape_adapted <- FALSE
  for (round in seq_len(1 + correlations)) {
    iter <- 1
    adapt_scale <- 1
    while ((round == 2 && !shape_adapted) ||
           (min(accRate) < min || max(accRate) > max) && iter <= max_iter) {
      message(date(), " Acceptance rate ", min(accRate),
              ", adapting ", ifelse(round == 1, "size", "shape"),
              " with scale ", adapt_scale)
      model <- output_to_proposal(adapt_wrapper, adapt_scale,
                                  correlations = (round == 2))
      add_options[["init-file"]] <- adapt_wrapper$result$output_file_name
      add_options[["init-np"]] <- samples - 1
      adapt_wrapper <-
        adapt_wrapper$clone(model = model, run = TRUE, add_options = add_options, ...)
      mcmc_obj <- mcmc(get_traces(adapt_wrapper))
      accRate <- acceptance_rate(adapt_wrapper)
      iter <- iter + 1
      if (min(accRate) < min) {
        adapt_scale <- adapt_scale / scale
      } else {
        adapt_scale <- adapt_scale * scale
      }
      if (round == 2) shape_adapted <- TRUE
     }
  }

  message(date(), " Acceptance rate:", min(accRate))

  if (iter > max_iter) {
    warning("Maximum of iterations reached")
  }

  return(adapt_wrapper)
}
