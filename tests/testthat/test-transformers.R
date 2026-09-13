test_that("trans_ratio scales by standard deviation", {
  x <- c(1, 2, 3, 4, 5)
  expect_equal(trans_ratio(x, "sd"), x / stats::sd(x))
})

test_that("trans_ratio scales by MAD", {
  x <- c(1, 2, 4, 8, 16)
  expect_equal(trans_ratio(x, "mad"), x / stats::mad(x))
})

test_that("trans_ratio validates the method", {
  expect_error(trans_ratio(1:5, "variance"))
})

test_that("trans_root applies the requested fractional power", {
  expect_equal(trans_root(c(1, 4, 9, 16)), c(1, 2, 3, 4))
  expect_equal(trans_root(c(1, 8, 27), fracRoot = 1 / 3), c(1, 2, 3), tolerance = 1e-12)
})

test_that("trans_log applies shift and base", {
  x <- c(0, 9, 99)
  expect_equal(trans_log(x, base = 10, shift = 1), c(0, 1, 2))
})
