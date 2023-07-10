context("Testing functions for adapting")

model <- system.file(package = "rbi", "PZ.bi")
example_output <- rbi::bi_read(
  system.file(package = "rbi", "example_output.nc")
)
bi <- rbi::libbi(model, nparticles = 1)
bi <- rbi::attach_data(bi, "output", example_output)
obs <- list(P_obs = data.frame(time = 0, value = 1.485121))

test_that("pMCMC adaptation works", {
  skip_on_cran()
  suppressWarnings(capture.output(
    prop_adapted <-
      adapt_proposal(bi,
        obs = obs, nsamples = 100, min = 0.5, max = 0.53,
        adapt = "both", seed = 3
      )
  ))
  prop_adapted2 <- adapt_proposal(
    prop_adapted, with = "transform-initial-to-param", min = 0.4, max = 0.5,
    adapt = "shape", seed = 5, scale = 0.5, quiet = TRUE
  )
  capture.output(part_adapted <- adapt_particles(prop_adapted, seed = 1))
  capture.output(part_adapted2 <- adapt_particles(bi, seed = 1, obs = obs))
  expect_true(part_adapted$options$nparticles > 0)
  expect_gt(length(bi_read(prop_adapted, file = "input")), 0)
  expect_gt(
    length(
      bi_read(adapt_proposal(prop_adapted, min = 0, max = 1), file = "input")
    ), 0
  )
})

test_that("warnings are given", {
  expect_warning(adapt_proposal(bi, size = TRUE, adapt = "shape"), "adapt size")
  expect_warning(
    adapt_proposal(bi, correlations = TRUE, adapt = "size"), "adapt shape"
  )
})

test_that("errors are reported", {
  expect_error(
    acceptance_rate(rbi::extract_sample(bi, 0), model = model),
    "just one sample"
  )
  expect_error(adapt_particles(bi, min = 1, max = 0), "greater")
  expect_error(adapt_proposal(bi, min = 1, max = 0), "max>min")
  expect_error(adapt_proposal(bi, adapt = "test"), "one of")
})

test_that("deprecated options are recognised", {
  expect_warning(adapt_proposal(bi, size = TRUE), "adapt")
  expect_warning(adapt_proposal(bi, correlations = TRUE), "adapt")
})
