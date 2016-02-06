##' Compute DIC for a libbi model
##'
##' @param read either a \code{libbi} object or a list of data frames, as returned by \code{bi_read}
##' @param burn number of iterations to discard as burn-in (if any)
##' @return DIC
##' @export
##' @author Sebastian Funk
compute_DIC <- function(read, burn)
{
    wrapper <- NULL
    
    ## read data
    if ("libbi" %in% class(read))
    {
        if (!read$run_flag)
        {
            stop("The model should be run first")
        }
        res <- bi_read(read)
        wrapper <- read
    } else if (is.data.frame(read))
    {
        res <- list(dummy = read)
    } else if (is.list(read))
    {
        res <- read
    } else
    {
        stop("'read' must be a 'libbi' object or a list of data frames or a data frame.")
    }

    ## convert to data.table
    res <- lapply(res, function(x) { if (is.data.frame(x)) { data.table(x) } else {x} })
    ## create a copy
    res <- lapply(res, copy)
    
    if (!missing(burn))
    {
        res <- lapply(res, function(x) {
            if ("np" %in% colnames(x)) { x <- x[np >= burn]}
        })
    }

    ## sample mean deviance
    mean_D <- -2 * mean(res[["loglikelihood"]]$value)

    ## effective number of parameters
    pd <- var(-2 * res[["loglikelihood"]]$value) / 2
    
    ## DIC
    return(mean_D + pd)
}

