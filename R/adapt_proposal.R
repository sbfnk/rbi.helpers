#' @rdname adapt_proposal
#' @name adapt_proposal
#' @title Adapt the proposal distribution of MCMC using the covariance
#'   of samples
#' @description This function takes the provided [rbi::libbi()] object and
#'   runs MCMC, adapting the proposal distribution until the desired
#'   acceptance rate is achieved. If a scale is given, it will be used
#'   to adapt the proposal at each iteration
#' @param x \code{link{libbi}} object
#' @param min minimum acceptance rate
#' @param max maximum acceptance rate
#' @param scale scale multiplier/divider for the proposal. If >1 this
#'   will be inverted.
#' @param max_iter maximum of iterations (default: 10)
#' @param adapt what to adapt; if "size" (default), the width of independent
#'   proposals will be adapted; if "shape", proposals will be dependent
#'   (following a multivariate normal) taking into account empirical
#'   correlations; if "both", the size will be adapted
#' before the shape
#' @param size (deprecated, use \code{{adapt}} instead) if TRUE (default:
#'   FALSE), the size of the (diagonal multivariate normal) proposal
#'   distribution will be adapted
#' @param correlations (deprecated, use \code{{adapt}} instead) if TRUE
#'   (default: FALSE), the shape of the (diagonal multivariate normal) proposal
#'   distribution will be adapted according to the empirical covariance
#' @param truncate if TRUE, the proposal distributions will be truncated
#'   according to the support of the prior distributions
#' @param quiet if set to TRUE, will not provide running output of particle
#'   numbers tested
#' @param ... parameters for [rbi::sample()]
#' @return a [rbi::libbi()] with the desired proposal distribution
#' @importFrom rbi get_traces sample enable_outputs get_dims
#'   attach_data
#' @examples
#' example_obs <- rbi::bi_read(system.file(package="rbi", "example_dataset.nc"))
#' example_model <- rbi::bi_model(system.file(package="rbi", "PZ.bi"))
#' example_bi <- rbi::libbi(model = example_model, obs = example_obs)
#' obs_states <- rbi::var_names(example_model, type="obs")
#' max_time <- max(vapply(example_obs[obs_states], function(x) {
#'   max(x[["time"]])
#' }, 0))
#' # adapt to acceptance rate between 0.1 and 0.5
#' \dontrun{
#'   adapted <- adapt_proposal(example_bi,
#'     nsamples = 100, end_time = max_time,
#'     min = 0.1, max = 0.5, nparticles = 256, correlations = TRUE
#'   )
#' }
#' @export
adapt_proposal <- function(x, min = 0, max = 1, scale = 2, max_iter = 10,
                           adapt = c("size", "shape", "both"), size = FALSE,
                           correlations = TRUE, truncate = TRUE, quiet = FALSE,
                           ...) {
  given_size <- NULL
  given_correlations <- NULL
  if (!missing(size)) {
    warning("'size' argument is deprecated; use 'adapt'")
    given_size <- size
  }
  if (!missing(correlations)) {
    warning("'correlations' argument is deprecated; use 'adapt'")
    given_correlations <- correlations
  }

  adapt <- match.arg(adapt)

  if (adapt == "both") adapt <- c("size", "shape")
  size <- ("size" %in% adapt)
  correlations <- ("shape" %in% adapt)
  if (!is.null(given_size) && size != given_size) {
    warning("'size' given but not in 'adapt'. Will not adapt size")
  }
  if (!is.null(given_correlations) && correlations != given_correlations) {
    warning(
      "'correlations' given but not in 'shape in 'adapt'. Will not adapt shape"
    )
  }

  if (min == 0 && max == 1) {
    return(x)
  }

  if (!(max > min)) stop("Must have max>min.")

  if (!quiet) message(date(), " Adapting the proposal distribution")

  ## ensure all parameters are saved to output file
  model_with_proposal <- update_proposal(
    x$model, truncate = truncate, blocks = c("parameter", "initial")
  )

  ## write initial covariance matrix
  adapted <- x
  adapted$model <- enable_outputs(model_with_proposal, type = "param")
  adaptation_vars <- get_mvn_params(adapted)
  adapted <- attach_data(
    adapted, file = "input", adaptation_vars, append = TRUE, overwrite = TRUE
  )

  if (!quiet) message(date(), " Initial trial run")
  adapted <- rbi::sample(adapted, ...)

  ## scale should be > 1 (it's a divider if acceptance rate is too
  ## small, multiplier if the acceptance Rate is too big)
  if (scale < 1) scale <- 1 / scale

  acc_rate <- acceptance_rate(adapted)
  rounds <- c()
  if (size) {
    rounds <- c(rounds, 1)
    shape_adapted <- FALSE
  } else {
    shape_adapted <- TRUE
  }
  if (correlations) rounds <- c(rounds, 2)
  for (round in rounds) {
    iter <- 1
    adapt_scale <- 1
    while ((round == 2 && !shape_adapted) ||
             (min(acc_rate) < min || max(acc_rate) > max ||
                !is.finite(acc_rate)) &&
               iter <= max_iter) {
      if (is.finite(acc_rate)) {
        if (!quiet) {
          message(
            date(), " Acceptance rate ", min(acc_rate),
            ", adapting ", ifelse(round == 1, "size", "shape"),
            " with scale ", adapt_scale
          )
        }
      }

      adaptation_vars <- get_mvn_params(
        adapted, correlations = (round == 2), scale = adapt_scale
      )
      adapted <-
        attach_data(adapted,
          file = "input", adaptation_vars, in_place = TRUE,
          overwrite = TRUE, quiet = TRUE
        )
      adapted <- rbi::sample(adapted, ...)
      acc_rate <- acceptance_rate(adapted)
      iter <- iter + 1
      if (min(acc_rate) < min) {
        adapt_scale <- adapt_scale / scale
      } else {
        adapt_scale <- adapt_scale * scale
      }
      if (round == 2) shape_adapted <- TRUE
      if (iter > max_iter) warning("Maximum of iterations reached")
    }
  }

  if (!quiet) message(date(), " Acceptance rate: ", min(acc_rate))

  ## put model back together (potential disabling output of some parameters)
  adapted$model <- model_with_proposal

  return(adapted)
}
