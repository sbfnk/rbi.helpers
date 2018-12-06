context("Testing functions for adapting")

model <- system.file(package="rbi", "PZ.bi")
example_output <- bi_read(system.file(package="rbi", "example_output.nc"))
bi <- libbi(model)
data <- bi_generate_dataset(bi, seed=1234)

test_that("the number of particles can be adapted",
{
  expect_true(adapt_particles(bi, obs=data, nsamples=100, seed=1234)$options$nparticles>0)
  expect_gt(length(bi_read(adapt_proposal(bi, obs=data, nsamples=100, min=0.1, max=0.52, adapt="both", seed=1234), file="input")), 0)
})
