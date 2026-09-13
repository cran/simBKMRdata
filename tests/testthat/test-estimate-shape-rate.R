test_that("estimate_mv_shape_rate computes method-of-moments parameters", {
  x <- data.frame(
    x1 = c(1, 2, 3, 4, 5),
    x2 = c(2, 4, 5, 8, 11)
  )
  out <- estimate_mv_shape_rate(x, using = "MoM")
  means <- colMeans(x)
  sds <- vapply(x, stats::sd, numeric(1))

  expect_equal(unname(out$shape_num), unname(means^2 / sds^2))
  expect_equal(unname(out$rate_num), unname(means / sds^2))
  expect_equal(out$sampSize, nrow(x))
  expect_equal(dim(out$sampCorr_mat), c(2, 2))
})

test_that("estimate_mv_shape_rate computes gMLE parameters for positive data", {
  x <- data.frame(
    x1 = c(1.2, 2.1, 3.7, 4.4),
    x2 = c(0.8, 1.5, 2.6, 5.1)
  )
  out <- estimate_mv_shape_rate(x, using = "gMLE")

  beta_scale <- colMeans(x * log(x)) - colMeans(x) * colMeans(log(x))
  alpha <- colMeans(x) / beta_scale
  expect_equal(unname(out$shape_num), unname(alpha))
  expect_equal(unname(out$rate_num), unname(1 / beta_scale))
})

test_that("estimate_mv_shape_rate validates inputs", {
  expect_error(estimate_mv_shape_rate(1:5), "data.frame or a matrix")
  expect_error(estimate_mv_shape_rate(data.frame(x = 1:4), using = "other"))
})
