test_that("generate_mvGamma_data returns positive data with expected dimensions", {
  set.seed(2026)
  out <- generate_mvGamma_data(
    sampSize = 40,
    sampCorr_mat = matrix(c(1, 0.4, 0.4, 1), 2, 2),
    shape_num = c(2, 3),
    rate_num = c(1, 2)
  )

  expect_s3_class(out, "data.frame")
  expect_equal(dim(out), c(40, 2))
  expect_true(all(is.finite(as.matrix(out))))
  expect_true(all(as.matrix(out) >= 0))
})

test_that("generate_mvGamma_data is reproducible under a fixed seed", {
  args <- list(
    sampSize = 10,
    sampCorr_mat = diag(2),
    shape_num = c(2, 2),
    rate_num = c(1, 1)
  )
  set.seed(11)
  a <- do.call(generate_mvGamma_data, args)
  set.seed(11)
  b <- do.call(generate_mvGamma_data, args)
  expect_equal(a, b)
})

test_that("generate_mvGamma_data validates parameter lengths", {
  expect_error(
    generate_mvGamma_data(10, diag(2), shape_num = 2, rate_num = c(1, 1)),
    "length of shape_num"
  )
  expect_error(
    generate_mvGamma_data(10, diag(2), shape_num = c(2, 2), rate_num = 1),
    "length of rate_num"
  )
})
