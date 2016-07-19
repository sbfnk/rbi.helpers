##' Merge parallel MCMC runs
##'
##' @param dir directory that contains the run files
##' @param pattern file pattern, excluding extensions
##' @param concatenate whether to concatenate traces, otherwise a list of traces (e.g., for \code{coda}) will be returned
##' @return traces
##' @author Sebastian Funk <sebastian.funk@lshtm.ac.uk>
merge_parallel_runs <- function(dir, pattern, concatenate = FALSE)
{
    mcmc_files <-
        sub("\\.rds$", "", list.files(dir, paste0(pattern, "\\.rds")))

    types <- c("obs", "init", "input")

    traces <- list()
    files <- list()
    model <- NULL

    for (file in paste(dir, mcmc_files, sep = "/"))
    {
        trace_file <- paste(file, "rds", sep = ".")
        if (file.exists(trace_file))
        {
            traces[[file]] <- readRDS(trace_file)
        }

        for (type in types)
        {
            filename <- paste0(file, "_", type, ".rds")
            if (is.null(files[[type]]) && file.exists(filename)) {
                files[[type]] <- readRDS(filename)
            }
        }

        model_filename <- paste(file, "bi", sep = ".")
        if (is.null(model) && file.exists(model_filename))
        {
            model <- bi_model(model_filename)
        }
    }

    if (concatenate & length(traces) > 0)
    {
        np_translate <- NULL
        combined <- lapply(names(traces[[1]]), function(x) {
            z <- rbindlist(lapply(seq_along(traces), function(y) {
                dt <- data.table(traces[[y]][[x]])
                dt[, unique_np := paste(np, y, sep = "_")]
            }))
            if (is.null(np_translate)) {
                np_translate <-
                    data.table(unique_np = unique(z$unique_np),
                               new_np = seq_along(unique(z$unique_np)) - 1)
            }
            z <- merge(z, np_translate, by = "unique_np")
            z[, np := NULL]
            z[, unique_np := NULL]
            setnames(z, "new_np", "np")
            setkey(z, np)
            z
        })

        if (length(combined) > 0) names(combined) <- names(traces[[1]])
        traces <- combined
    }

    ret <- list()

    if (length(traces) > 0) ret <- c(ret, list(traces = traces))
    if (length(files) > 0) ret <- c(ret, files)
    if (!is.null(model)) ret <- c(ret, list(model = model))

    return(ret)
}
