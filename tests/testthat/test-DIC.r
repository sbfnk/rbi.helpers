context("Testing functions for information criteria")

example_run <- rbi::bi_read(system.file(package = "rbi", "example_output.nc"))
example_model_file <- system.file(package = "rbi", "PZ.bi")
example_bi <- rbi::attach_data(
  rbi::libbi(example_model_file), "output", example_run
)

test_that("We can calculate DIC", {
  expect_equal(round(DIC(example_bi, burn = 10)), 174)
})
