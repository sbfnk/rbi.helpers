#' @rdname sample_observations
#' @name sample_observations
#' @title Sample observations from trajectories of a libbi model, one at each time step
#' @param wrapper a \code{\link{libbi}} object that has been run
#' @param read_options options for \code{\link{bi_read}}
#' @param ... further options to be passed to \code{libbi} (with run = TRUE)
#' @return \code{\link{libbi}} object with sampled observations (if an \code{output_file_name} has been given), or a list of data frames containing the results
#' @export
#' 
sample_observations <- function(wrapper, read_options, ...){
  if (!wrapper$run_flag) {
    stop("The model should be run first")
  }

  times <- bi_read(wrapper, "time", vector = TRUE)
  iterations <- bi_dim_len(wrapper$result$output_file_name, "np")
  end_time <- max(times) + 1
  start_time <- min(times)

  model <- wrapper$model$clone()
  model$remove_block("transition")
  model$remove_block("parameter")
  model$remove_block("initial")

  run_joint <- wrapper$clone(run = TRUE, model = model,
                             "start-time" = start_time,
                             "end-time" = end_time,
                             "noutputs" = end_time - start_time,
                             target = "joint",
                             nsamples = iterations,
                             input = wrapper, init = wrapper,
                             ...)

  if (missing(read_options)) read_options <- NULL
  res <- do.call(bi_read, c(list(read = run_joint,
                                 vars = model$get_vars("obs")),
                            read_options))
  if (is.data.frame(res)) {
    res <- list(res)
    names(res) <- model$get_vars("obs")
  }
  for (var in names(res)) {
      if ("nr" %in% names(res[[var]])) {
          res[[var]]$nr <- res[[var]]$nr - 1
          res[[var]] <- res[[var]][res[[var]]$nr >= 0, ]
      }
  }
  ## if an output file has been given, we re-write
  args <- match.call()
  if ("output_file_name" %in% names(args))
  {
      bi_write(run_joint$result$output_file_name, res)
      return(run_joint)
  } else
  {
      return(res)
  }
}
