#' @rdname acceptance_rate
#' @name acceptance_rate
#' @title Compute acceptance rate
#' @description
#' This function takes the provided \code{\link{libbi}} object which has been
#' run, or a bi file, and returns a the acceptance rate
#' @param ... parameters to \code{\link{get_traces}} (especially 'run')
#' @return acceptance rate
#' @importFrom coda mcmc rejectionRate
#' @importFrom rbi get_traces
#' @export
#' @examples
#' example_run_file <- system.file(package="rbi.helpers", "example_run.nc")
#' example_model_file <- system.file(package="rbi", "PZ.bi")
#' example_bi <- add_output(libbi(example_model_file), example_run_file)
#' acceptance_rate(example_bi)
#' @author Sebastian Funk
acceptance_rate <- function(...) {

  mcmc_obj <- coda::mcmc(rbi::get_traces(...))
  if (nrow(mcmc_obj) > 1) {
    accRate <- max(1 - coda::rejectionRate(mcmc_obj))
  } else {
    stop("Cannot compute acceptance rate from just one sample. Try setting 'nsamples'")
  }

  return(accRate)
}
