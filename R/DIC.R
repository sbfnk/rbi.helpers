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
##' @param burn number of iterations to discard as burn-in (if any)
##' @param bootstrap number of bootstrap samples to take, 0 to just take data
##' @param ... not used
##' @return DIC
##' @import data.table
##' @importFrom stats var
##' @export
##' @examples
##' example_run <- bi_read(system.file(package="rbi", "example_output.nc"))
##' example_model_file <- system.file(package="rbi", "PZ.bi")
##' example_bi <- attach_data(libbi(example_model_file), "output", example_run)
##' DIC(example_bi)
##' @author Sebastian Funk
DIC.libbi <- function(x, burn, bootstrap = 0, ...)
{
    res <- bi_read(x)

    ## convert to data.table
    res <- lapply(res, function(x) { if (is.data.frame(x)) { data.table::data.table(x) } else {x} })

    if (!missing(burn))
    {
        burned <- lapply(res, function(x) {
            if ("np" %in% colnames(x)) {
                x <- x[get("np") >= burn]
            }
        })
    } else
    {
        burned <- res
    }

    ## sample mean deviance
    mean_D <- mean(-2 * burned[["loglikelihood"]]$value)

    ## effective number of parameters
    pd <- stats::var(-2 * burned[["loglikelihood"]]$value) / 2

    ## DIC
    return(mean_D + pd)
}

