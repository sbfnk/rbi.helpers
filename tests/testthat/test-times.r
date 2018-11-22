context("Testing functions for times operations")

example_run <- bi_read(system.file(package="rbi", "example_output.nc"))
example_model_file <- system.file(package="rbi", "PZ.bi")
example_bi <- attach_data(libbi(example_model_file), "output", example_run)

test_that("We can convert numeric to actual times",
{
  expect_equal(class(numeric_to_time(example_bi, origin=as.Date("2018-01-01"), unit="day")$Z$time), "Date")
  expect_true("POSIXt" %in%  class(numeric_to_time(bi_read(example_bi), origin=as.Date("2018-01-01"), unit="2 minutes")$Z$time))
})

test_that("We can convert actual to numeric times",
{
  out <- numeric_to_time(example_bi, origin=as.Date("2018-01-01"), unit="day")
  expect_equal(class(time_to_numeric(out, origin=as.Date("2018-01-01"), unit="day")$Z$time), "numeric")
  expect_equal(class(time_to_numeric(out$Z, origin=as.Date("2018-01-01"), unit="day")$time), "numeric")
})

test_that("Errors are recognised",
{
  expect_error(time_to_numeric(3), "must be")
  expect_error(numeric_to_time(3), "must be")
})
