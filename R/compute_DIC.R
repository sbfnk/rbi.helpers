##' Compute DIC for a libbi model
##'
##' @param read either a \code{libbi} object or a list of data frames, as returned by \code{bi_read}
##' @param model model file or a \code{bi_model} object (if \code{read} is not a \code{libbi} object)
##' @param burn number of iterations to discard as burn-in (if any)
##' @param ... options for \link{\code{libbi}}
##' @return DIC
##' @export
##' @author Sebastian Funk
compute_DIC <- function(read, model, burn, ...)
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
        res <- lapply(res, function(x) { if ("np" %in% colnames(x)) { x <- x[np >= burn]}  })
    }

    ## read model
    if (missing(model))
    {
        if ("libbi" %in% class(read))
        {
            model <- read$model
        } else
        {
            stop("If 'read' is not a libbi object, 'model' must be given.")
        }
    } else
    {
        if ("libbi" %in% class(read))
        {
            warning("'model' overwrites the model given in 'read'.x")
        }
        if (is.character(model)) {
            model <- bi_model(model)
        } else if (!("bi_model" %in% class(model))) {
            stop("'model' must be either a 'bi_model' object or a path to a valid model file in LibBi's syntax")
        }
    }

    ## read parameters
    parameters <- model$get_vars("param")

    ## work out mean parameter
    theta_mean <- list()
    for (param in parameters)
    {
        theta_mean[[param]] <- mean(res[[param]][, value])
    }

    ## set up libbi to work out likelihood at mean parameter
    if (is.null(wrapper))
    {
        wrapper <- libbi(client = "sample", model = model, init = theta_mean, ...)
    } else
    {
        if ("init-file" %in% names(wrapper$global_options))
        {
            init <- bi_read(wrapper$global_options[["init-file"]])
            for (param in parameters)
            {
                init[[param]] <- theta_mean[[param]]
            }
        }
        if ("nsamples" %in% names(wrapper$global_options))
        {
            global_options[["nsamples"]] <- 1
        }
    }

    wrapper$run(target = "posterior")

    ## likelihood at mean parameter
    ll_mean <- bi_read(wrapper)$loglikelihood

    ## deviance at mean parameter
    D_mean <- -2 * ll_mean

    ## sample mean deviance
    mean_D <- -2 * mean(res[["loglikelihood"]]$value)

    ## effective number of parameters
    pd <- mean_D - D_mean
    
    ## DIC
    return(D_mean + 2 * pd)
}

