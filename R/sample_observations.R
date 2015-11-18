#' @rdname sample_observations
#' @name sample_observations
#' @title Sample observations from trajectories of a libbi model, one at each time step
#' @param wrapper a \code{\link{libbi}} object that has been run
#' @param output_file_name the name of the output file to write
#' @param ... further options to be passed to \code{libbi} (with run = TRUE)
#' @export
#' 
sample_observations <- function(wrapper, output_file_name, ...){
  if (!wrapper$run_flag) {
    stop("The model should be run first")
  }

  times <- bi_read(wrapper$result$output_file_name, "time", vector = TRUE)
  iterations <- bi_dim_len(wrapper$result$output_file_name, "np")
  end_time <- max(times) + 1
  start_time <- min(times)

  model <- wrapper$model$clone()
  model$remove_block("transition")
  model$remove_block("parameter")
  model$remove_block("initial")

  run_joint <- wrapper$clone(run = TRUE, model = model, 
                             add_options = list("start-time" = start_time, 
                                                "end-time" = end_time, 
                                                "noutputs" = end_time - start_time, 
                                                target = "joint",
                                                nsamples = iterations),
                             input = wrapper, init = wrapper,
                             output_file_name = output_file_name, ...)

  res <- bi_read(output_file_name)
  for (var in names(res)) {
      if ("nr" %in% names(res[[var]])) {
          res[[var]]$nr <- res[[var]]$nr - 1
          res[[var]] <- res[[var]][res[[var]]$nr >= 0, ]
      }
  }
  bi_write(output_file_name, res)
}
