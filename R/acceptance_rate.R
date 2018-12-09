#' @rdname acceptance_rate
#' @name acceptance_rate
#' @title Compute acceptance rate
#' @description
#' This function takes the provided \code{\link{libbi}} object which has been
#' run, or a bi file, and returns a the acceptance rate
#' @param ... parameters to \code{\link{get_traces}} (especially 'x')
#' @return acceptance rate
#' @importFrom rbi get_traces
#' @export
#' @examples
#' example_run <- bi_read(system.file(package="rbi.helpers", "example_run.nc"))
#' example_model_file <- system.file(package="rbi", "PZ.bi")
#' example_bi <- attach_data(libbi(example_model_file), "output", example_run)
#' acceptance_rate(example_bi)
#' @author Sebastian Funk
acceptance_rate <- function(...) {

  ## we use the same method as the rejectionRate function of the`coda` library (v. 0.10-7):
  ## Martyn Plummer, Nicky Best, Kate Cowles and Karen Vines (2006). CODA:
  ## Convergence Diagnosis and Output Analysis for MCMC, R News, vol 6,
  ## 7-11

  x <- rbi::get_traces(...)
  if (nrow(x) > 1) {
    accRate <- 1 - mean(apply(x[-nrow(x), , drop = FALSE] == x[-1,, drop = FALSE], 1, all))
  } else {
    stop("Cannot compute acceptance rate from just one sample. Try setting 'nsamples'")
  }

  return(accRate)
}
