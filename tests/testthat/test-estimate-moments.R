test_that("estimate_mv_moments returns expected summaries", {
  x <- data.frame(
    x1 = c(1, 2, 3, 4),
    x2 = c(2, 4, 6, 8)
  )
  out <- estimate_mv_moments(x)

  expect_named(out, c("sampSize", "mean_vec", "sampSD", "sampCorr_mat", "sampSkew"))
  expect_equal(out$sampSize, 4)
  expect_equal(unname(out$mean_vec), c(2.5, 5))
  expect_equal(unname(out$sampSD), c(stats::sd(x$x1), stats::sd(x$x2)))
  expect_equal(unname(out$sampCorr_mat), matrix(c(1, 1, 1, 1), 2, 2))
  expect_length(out$sampSkew, 2)
})

test_that("estimate_mv_moments accepts matrices", {
  x <- matrix(c(1, 2, 3, 3, 2, 1), ncol = 2)
  expect_equal(
    estimate_mv_moments(x)$mean_vec,
    estimate_mv_moments(as.data.frame(x))$mean_vec
  )
})

test_that("estimate_mv_moments handles missing values in basic moments", {
  x <- data.frame(x1 = c(1, 2, NA, 4), x2 = c(2, 3, 4, 5))
  out <- suppressWarnings(estimate_mv_moments(x))
  expect_equal(unname(out$mean_vec), c(mean(x$x1, na.rm = TRUE), mean(x$x2)))
  expect_equal(unname(out$sampSD), c(stats::sd(x$x1, na.rm = TRUE), stats::sd(x$x2)))
})

test_that("estimate_mv_moments returns NA skewness for constant columns", {
  x <- data.frame(x1 = rep(2, 5), x2 = 1:5)
  out <- suppressWarnings(estimate_mv_moments(x))
  expect_true(is.na(out$sampSkew[1]))
})

test_that("estimate_mv_moments rejects unsupported inputs", {
  expect_error(estimate_mv_moments(1:5), "data.frame or a matrix")
})
