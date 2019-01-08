#' @rdname get_mvn_params
#' @name get_mvn_params
#' @title Construct a covariance matrix
#' @description
#' This function takes the provided \code{\link{libbi}} which has been
#' run and returns the square root of the covariance matrix, which
#' can be used for proposal distributions
#' @param x a \code{\link{libbi}} which has been run
#' @param scale a factor by which to scale all elements of the covariance matrix
#' @param correlations logical; if TRUE, correlations are taken into account when constructing the parameters
#' @param fix if this is set, all elements of the covariance matrix will be set to it
#' @inheritParams update_proposal
#' @importFrom data.table setnames data.table
#' @importFrom reshape2 melt
#' @importFrom stats cov sd
#' @importFrom Matrix nearPD
#' @return the updated bi model
#' @keywords internal
get_mvn_params <- function(x, scale=1, correlations=TRUE, fix) {

  ## get parameters from update_proposal
  blocks <- c("parameter", "initial")
  params <- lapply(blocks, function(block) {
    block_lines <-
      grep("~", get_block(x$model, paste("proposal", block, sep="_")), value=TRUE)
    sub("[[:space:]]*~.*$", "", block_lines)
  })
  names(params) <- blocks
  dimless_params <- lapply(params, function(param) {
    sub("[[:space:]]*\\[.*$", "", param)
  })
  names(dimless_params) <- names(params)

  dims <- get_dims(x$model)

  vars <- list()

  for (block in names(params)) {
    ## read parameters

    if (!missing(fix)) {
      A <- matrix(fix, ncol=length(params[[block]]), nrow=length(params[[block]]))
      rownames(A) <- params[[block]]
      colnames(A) <- params[[block]]
    } else if (x$run_flag){
      res <- bi_read(x, vars=unique(dimless_params[[block]]), init_to_param = TRUE)

      l <- list()
      for (param in names(res))
      {
        if (prod(dim(res[[param]])) > 0)
        {
          y <- data.table(res[[param]])
          ## extract columns that are dimensions
          dim_colnames <- setdiff(colnames(y), c("np", "value"))
          unique_dims <- unique(y[, dim_colnames, with=FALSE])
          if (sum(dim(unique_dims)) > 0)
          {
            ## for parameters with dimensions, create a parameter for each
            ## possible dimension(s)
            a <- lapply(seq_len(nrow(unique_dims)), function(row) {
              merge(unique_dims[row, ], y)
            })
            ## convert factor to integer
            for (dim_colname in dim_colnames) {
              if (is.factor(unique_dims[[dim_colname]])) {
                unique_dims[[dim_colname]] <-
                  as.integer(unique_dims[[dim_colname]]) - 1
              }
            }
            ## create correct parameter names (including the dimensions)
            if (length(a) > 0)
            {
              names(a) <- unname(apply(unique_dims, 1, function(row) {
                paste0(param, "[", paste(row, collapse = ","), "]")
              }))
            }
          } else
          {
            a <- list(y)
            names(a) <- param
          }
        }

        ## loop over all parameters (if dimensionsless, just the parameter,
        ## otherwise all possible dimensions) and remove all dimension columns
        a <- lapply(names(a), function(p) {
          for (col in colnames(unique_dims)) {
            a[[p]][[col]] <- NULL
          }
          data.table::setnames(a[[p]], "value", p)
        })
        l <- c(l, a)
      }

      if (length(l) > 0) {
        ## create a wide table of all the parameters, for calculating
        ## covariances
        wide <- l[[1]]
        if (length(l) > 1) {
          for (i in seq(2, length(l))) {
            wide <- merge(wide, l[[i]], by=intersect(colnames(wide), colnames(l[[i]])))
          }
        }
        wide[["np"]] <- NULL

        if (ncol(wide) > 0)
        {
          if (ncol(wide) > 1 && correlations) {
            ## calculate the covariance matrix
            c <- 2.38**2 * stats::cov(wide) / ncol(wide)
            A <- tryCatch(t(as.matrix(chol(nearPD(c)$mat))), error=function(e) NULL)
          } else {
            A <- NULL
          }

          if (is.null(A)) {
            A <- diag(apply(wide, 2, sd))
            rownames(A) <- colnames(wide)
            colnames(A) <- colnames(wide)
          }

          dA <- diag(A)
          for (i in seq_along(dA)) {
            if (dA[i]==0) {
              A[i, i] <- abs(mean(wide[, i, with=FALSE][[1]]))/10
            }
          }
        } else  {
          A <- matrix(ncol=0, nrow=0)
        }
      } else {
        A <- NULL
      }
    } else {
      A <- NULL
    }

    if (is.null(A)) {
      A <- diag(ncol=length(params[[block]]), nrow=length(params[[block]]))
      rownames(A) <- params[[block]]
      colnames(A) <- params[[block]]
    }

    if (prod(dim(A) > 0)) {
      if (!missing(scale)) {
        A <- A * sqrt(scale)
      }
      ## generate variables
      cov_input_name <- paste("__proposal", block, "cov", sep="_")
      dim_name <- paste("__dim", block, "cov", sep="_")
      long_cov <- reshape2::melt(A, varnames=paste(dim_name, 1:2, sep="."))
      vars[[cov_input_name]] <- long_cov
    }
  }

  return(vars)
}
