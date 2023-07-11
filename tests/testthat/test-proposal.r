context("Testing generating multivariate proposals")

set.seed(12124)

model_str <- "
model test {
  const no_a = 2
  const no_b = 3
  dim a(no_a)
  dim b(no_b)

  obs M[a]

  state N[a] (has_input = 0)
  noise e[a, b]
  param m[a, b]
  param p

  sub parameter {
    m[0,0] ~ beta()
    m[0,1] ~ gamma()
    m[0,2] ~ gaussian(mean = 1)
    m[1,b] ~ truncated_gaussian(mean = 1, lower = m[0,b])
    p ~ uniform(lower = 0, upper = 1)
  }

  sub initial {
    N[a] <- 1
  }

  sub transition {
    e[a, b] ~ gaussian(mean = m[a,b])
    N[a] <- N[a] + e[a, 0] +
       e[a, 1]
  }

  sub observation {
    inline x = m
    M[a] ~ gaussian(mean = N[a])
  }
}
"

model <- rbi::bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi <- rbi::libbi(model, dims = list(a = c("first", "second")))
test_output <- Map(
  function(x) {
    if (is.data.frame(x)) x$value <- abs(rnorm(nrow(x)))
    x
  },
  list(
    e = data.frame(
      expand.grid(time = 0:1, a = c("first", "second"), b = 0:2, np = 0:1)
    ),
    N = data.frame(
      expand.grid(time = 0:1, a = c("first", "second"), np = 0:1)
    ),
    m = data.frame(expand.grid(a = c("first", "second"), b = 0:2, np = 0:1)),
    p = data.frame(np = 0:1),
    close = 0,
    loglikelihood = data.frame(np = 0:1),
    logprior = data.frame(np = 0:1)
  )
)
bi <- rbi::attach_data(bi, "output", test_output)
test_that("multivariate proposals can be generated with multiple parameters", {
  skip_on_cran()

  updated_model <- update_proposal(model)
  bi$model <- updated_model

  expect_gt(length(get_mvn_params(bi)), 0)
  expect_gt(length(get_mvn_params(bi, fix = 0)), 0)
})

test_that("multivariate proposals can be generated with one parameter", {

  skip_on_cran()

  updated_model <- update_proposal(model)
  bi$model <- updated_model

  bi$model <- rbi::fix(bi$model, m = 0)
  test_output$m <- NULL
  test_output$p$value <- 1 ## no variation
  bi <- rbi::attach_data(bi, "output", test_output)
  expect_gt(length(get_mvn_params(bi)), 0)

  test_output$p$value <- 1 ## no variation
  expect_gt(length(get_mvn_params(bi)), 0)
})

test_that("multivariate proposals can be generated if covariance is 0", {

  skip_on_cran()

  updated_model <- update_proposal(model)
  bi$model <- updated_model

  test_output$m$value <- 1
  test_output$p$value <- 1
  bi <- rbi::attach_data(bi, "output", test_output)

  expect_gt(length(get_mvn_params(bi)), 0)
})

test_that("errors are reported", {
  expect_error(update_proposal(3), "must be a")
})
