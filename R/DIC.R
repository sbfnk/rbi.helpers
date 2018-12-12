#' @export
DIC <- function(x, ...) UseMethod("DIC")
##' @name DIC
##' @rdname DIC
##' @title Compute Deviance Information Criterion (DIC) for a libbi model
##'
##' @description
##' Computes the DIC of a libbi object containing Monte-Carlo samples. The effective number of parameters is calculated following Gelman et al., Bayesian Data Analysis: Second Edition, 2004, p. 182.
##'
##' @param x a \code{libbi} object
##' @param bootstrap number of bootstrap samples to take, 0 to just take data
##' @param ... any parameters to be passed to `bi_read` (e.g., `burn`)
##' @return DIC
##' @importFrom stats var
##' @export
##' @examples
##' example_run <- bi_read(system.file(package="rbi", "example_output.nc"))
##' example_model_file <- system.file(package="rbi", "PZ.bi")
##' example_bi <- attach_data(libbi(example_model_file), "output", example_run)
##' DIC(example_bi)
##' @author Sebastian Funk
DIC.libbi <- function(x, bootstrap = 0, ...)
{
    res <- bi_read(x, ...)

    ## sample mean deviance
    mean_D <- mean(-2 * res[["loglikelihood"]]$value)

    ## effective number of parameters
    pd <- stats::var(-2 * res[["loglikelihood"]]$value) / 2

    ## DIC
    return(mean_D + pd)
}

