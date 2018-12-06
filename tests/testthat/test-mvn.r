context("Testing generating multivariate proposals")

model_str <- "
model test {
  const no_a = 2
  const no_b = 2
  dim a(no_a)
  dim b(no_b)

  obs M[a]

  state N[a] (has_input = 0)
  noise e[a, b]
  param m[a, b]

  sub parameter {
    m[a,b] ~ truncated_gaussian(lower=0)
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

model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi <- libbi(model, dims=list(a=c("first", "second")))
test_output <-
    Map(
        function(x) { if (is.data.frame(x)) x$value <- abs(rnorm(nrow(x))); x },
        list(e=data.frame(expand.grid(time=0:1, a=c("first", "second"), b=0:1, np=0:1)),
             N=data.frame(expand.grid(time=0:1, a=c("first", "second"), np=0:1)),
             m=data.frame(expand.grid(a=c("first", "second"), b=0:1, np=0:1)),
             close=0,
             loglikelihood=data.frame(np=0:1),
             logprior=data.frame(np=0:1)))
bi <- attach_data(bi, "output", test_output)

test_that("multivariate proposals can be generated", 
{
  updated_model <- update_proposal(model)
  bi$model <- updated_model

  expect_gt(length(get_mvn_params(bi)), 0)
  expect_gt(length(get_mvn_params(bi, fix=0)), 0)

  test_output2 <- test_output
  test_output2$m$value <- 1
  bi2 <- attach_data(bi, "output", test_output2)

  expect_gt(length(get_mvn_params(bi2)), 0)

})

test_that("errors are reported",
{
    expect_error(update_proposal(3), "must be a")
})

test_that("deprecated options are recognised",
{
  expect_warning(adapt_proposal(bi, size=TRUE), "adapt")
  expect_warning(adapt_proposal(bi, correlations=TRUE), "adapt")
})
