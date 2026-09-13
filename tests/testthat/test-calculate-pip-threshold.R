test_that("calculate_pip_threshold reproduces the documented value", {
  out <- calculate_pip_threshold(absCV = 7.5, sampSize = 300)
  expect_equal(out, 0.6549943, tolerance = 1e-7)
})

test_that("calculate_pip_threshold can calculate inputs from y", {
  y <- c(2, 3, 4, 5, 6, 7)
  explicit <- calculate_pip_threshold(
    absCV = abs(stats::sd(y) / mean(y)),
    sampSize = length(y)
  )
  from_y <- calculate_pip_threshold(y = y)
  expect_equal(from_y, explicit)
})

test_that("calculate_pip_threshold requires both summary inputs", {
  expect_error(
    calculate_pip_threshold(absCV = 1),
    "both absCV and sampSize are required"
  )
  expect_error(
    calculate_pip_threshold(sampSize = 100),
    "both absCV and sampSize are required"
  )
})

test_that("calculate_pip_threshold accepts user-supplied coefficients", {
  coeffs <- list(A = 0, K = 1, C = 1, betaAbsCV = 0, betaSampSize = 1)
  out <- calculate_pip_threshold(absCV = 2, sampSize = 10, coeffs_ls = coeffs)
  expect_equal(out, 0.5)
})
