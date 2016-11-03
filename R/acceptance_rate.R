#' @rdname acceptance_rate
#' @name acceptance_rate
#' @title Compute acceptance rate
#' @description
#' This function takes the provided \code{\link{libbi}} object which has been
#' run, or a bi file, and returns a the acceptance rate
#' @param ... parameters to \code{\link{get_traces}} (especially 'run')
#' @importFrom coda mcmc rejectionRate
#' @return acceptance rate
#' @export
acceptance_rate <- function( ...) {

  mcmc_obj <- coda::mcmc(get_traces(...))
  accRate <- max(1 - coda::rejectionRate(mcmc_obj))

  return(accRate)
}
