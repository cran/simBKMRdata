test_that("simulate_group_gamma combines groups", {
  params <- list(
    Male = list(sampSize = 6, sampCorr_mat = diag(2), shape_num = c(2, 2), rate_num = c(1, 1)),
    Female = list(sampSize = 4, sampCorr_mat = diag(2), shape_num = c(3, 4), rate_num = c(2, 2))
  )
  set.seed(12)
  out <- simulate_group_gamma(params, "Sex")

  expect_equal(nrow(out), 10)
  expect_equal(ncol(out), 3)
  expect_equal(as.integer(table(out$Sex)), c(4L, 6L))
})

test_that("simulate_group_gamma checks names and required parameters", {
  expect_error(
    simulate_group_gamma(list(list(sampSize = 2, sampCorr_mat = diag(2), shape_num = c(2, 2), rate_num = c(1, 1))), "grp"),
    "named sublists"
  )
  bad <- list(A = list(sampSize = 2, sampCorr_mat = diag(2), shape_num = c(2, 2)))
  expect_error(simulate_group_gamma(bad, "grp"), "must contain")
})

test_that("simulate_group_gaussian combines groups", {
  params <- list(
    Male = list(sampSize = 6, mean_vec = c(0, 1), sampSD = c(4, 7), sampCorr_mat = diag(2)),
    Female = list(sampSize = 4, mean_vec = c(2, 3), sampSD = c(8, 9), sampCorr_mat = diag(2))
  )
  set.seed(13)
  out <- simulate_group_gaussian(params, "Sex")

  expect_equal(nrow(out), 10)
  expect_equal(ncol(out), 3)
  expect_equal(as.integer(table(out$Sex)), c(4L, 6L))
})

test_that("simulate_group_gaussian retains the documented unit-variance behavior", {
  params <- list(
    A = list(
      sampSize = 5000,
      mean_vec = c(0, 0),
      sampSD = c(10, 20),
      sampCorr_mat = matrix(c(1, 0.35, 0.35, 1), 2, 2)
    )
  )
  set.seed(2026)
  out <- simulate_group_gaussian(params, "grp")

  # sampCorr_mat is intentionally passed directly to MASS::mvrnorm() as Sigma.
  expect_equal(stats::sd(out$V1), 1, tolerance = 0.06)
  expect_equal(stats::sd(out$V2), 1, tolerance = 0.06)
  expect_equal(stats::cor(out$V1, out$V2), 0.35, tolerance = 0.05)
})

test_that("simulate_group_gaussian checks names and required parameters", {
  expect_error(
    simulate_group_gaussian(list(list(sampSize = 2, mean_vec = c(0, 0), sampSD = c(1, 1), sampCorr_mat = diag(2))), "grp"),
    "named sublists"
  )
  bad <- list(A = list(sampSize = 2, mean_vec = c(0, 0), sampCorr_mat = diag(2)))
  expect_error(simulate_group_gaussian(bad, "grp"), "must contain")
})
