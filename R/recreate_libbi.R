#' @rdname recreate_libbi
#' @name recreate_libbi
#' @title Recreate libbi wrapper
#' @description
#' This function links a model and output file name of a previous run and
#' returns a \code{\link{libbi}} object, so as to be able to pass it to analysis
#' functions
#'
#' @param model either a \code{\link{bi_model}} object, or a path to a model file
#' @param output_file_name the NetCDF file holding the output of the \code{libbi} run
#' @param client defaults to "sample"; this won't be needed in most cases but prevents a warning.
#' @param ... further options to pass to \code{libbi}, e.g. global_options
#' @return \code{\link{libbi}} object that can be analysed
#' @importFrom rbi libbi
#' @examples
#' example_run_file <- system.file(package="rbi.helpers", "example_run.nc")
#' example_model_file <- system.file(package="rbi", "PZ.bi")
#' example_bi <- recreate_libbi(example_model_file, example_run_file)
#' @export
recreate_libbi <- function(model, output_file_name, client = "sample", ...)
{

  wrapper <- rbi::libbi(client = "sample", model = model, ...)
  wrapper$output_file_name <- output_file_name
  wrapper$run_flag <- TRUE

  return(wrapper)
}
