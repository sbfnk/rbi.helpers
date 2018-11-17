#' @rdname output_to_cov
#' @name output_to_cov
#' @title Construct a covariance matrix
#' @description
#' This function takes the provided \code{\link{libbi}} which has been
#' run and returns the square root of the covariance matrix, which
#' can be used for proposal distributions
#' @param x a \code{\link{libbi}} which has been run
#' @param scale a factor by which to scale all the standard deviations
#' @inheritParams update_proposal
#' @importFrom data.table setnames melt data.table
#' @importFrom stats cov
#' @importFrom Matrix nearPD
#' @return the updated bi model
#' @keywords internal
output_to_cov <- function(x, scale, correlations=FALSE) {

  if (!x$run_flag) {
    stop("The model should be run first")
  }

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

  ## read parameters
  res <- bi_read(x, vars=unlist(unique(dimless_params)), init.to.param = TRUE)

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

  ## create a wide table of all the parameters, for calculating
  ## covariances
  wide <- l[[1]]
  if (length(l) > 1) {
    for (i in seq(2, length(l))) {
      wide <- merge(wide, l[[i]], by=intersect(colnames(wide), colnames(l[[i]])))
    }
  }
  wide[["np"]] <- NULL

  vars <- list()

  for (block in names(params))
  {
    wide_block <- wide[, params[[block]], with=FALSE]

    if (ncol(wide_block) > 0)
    {
      if (ncol(wide_block) > 1 && correlations) {
        ## calculate the covariance matrix
        c <- 2.38**2 * stats::cov(wide_block) / ncol(wide_block)
        tryCatch(A <- t(as.matrix(chol(nearPD(c)$mat))), error=function(e) NULL)
        if (is.null(A))
        {
          warning("Cholesky decomposition failed; ",
                  "will try to adapt via independent univariate sampling first.")
        }
      } else {
        A <- NULL
      }

      if (is.null(A)) {
        A <- diag(apply(wide_block, 2, sd))
        rownames(A) <- params[[block]]
        colnames(A) <- params[[block]]
      }

      if (!missing(scale)) {
        A <- A * scale
      }

      for (i in seq_len(ncol(A))) {
        if (A[i, i]==0) {
          A[i, i] <- mean(wide_block[, i, with=FALSE][[1]])/10
        }
      }

      ## generate variables - first check if covariate input exists in data
      cov_input_name <- paste("__proposal", block, "cov", sep="_")
      input_var_names <- var_names(x$model, type="input")
      if (cov_input_name %in% input_var_names) {
        dim_name <- paste("__dim", block, "cov", sep="_")
        long_cov <- data.table::melt(A, varnames=paste(dim_name, 1:2, sep="."))
        vars[[cov_input_name]] <- long_cov
      } else if (correlations) {
        stop("'correlations=TRUE' but model file doesn't ask for correlation input.")
      } else {
        unique_params <- unique(dimless_params[[block]])
        std_vars <- lapply(unique_params, function(param)
        {
          dims <- get_dims(x$model)
          var_dims <-
            sub("^.*\\[(.*)\\]", "\\1", var_names(x$model, vars=param, dim=TRUE))
          if (var_dims == param) {
            ## no dimensions
            A[param, param]
          } else {
            var_dim_names <- split(var_dims, ",")[[1]]
            ids <- which(dimless_params[[block]] == param)
            std <- diag(A[ids, ids])
            std_dims <- sub("^.*\\[(.*)\\]", "\\1", names(std))
            df <-
              as.data.frame(do.call(rbind, lapply(strsplit(std_dims, ","), as.integer)))
            colnames(df) <- var_dim_names
            cbind(df, data.frame(value=unname(std)))
          }
        })
        names(std_vars) <- unique_params
        for (name in names(std_vars)) {
          std_input_name <- paste("__std", name, sep="_")
          vars[[std_input_name]] <- std_vars[[name]]
        }
      }
    }
  }

  return(vars)
}
