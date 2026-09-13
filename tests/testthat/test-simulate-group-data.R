test_that("simulate_group_data generates grouped Gamma data", {
  params <- list(
    A = list(sampSize = 5, sampCorr_mat = diag(2), shape_num = c(2, 3), rate_num = c(1, 2)),
    B = list(sampSize = 7, sampCorr_mat = diag(2), shape_num = c(3, 4), rate_num = c(2, 2))
  )
  set.seed(1)
  out <- simulate_group_data(params, generate_mvGamma_data, "grp")

  expect_equal(nrow(out), 12)
  expect_equal(ncol(out), 3)
  expect_equal(as.integer(table(out$grp)), c(5L, 7L))
  expect_identical(rownames(out), as.character(seq_len(nrow(out))))
})

test_that("simulate_group_data generates grouped standardized Gaussian data", {
  params <- list(
    A = list(sampSize = 8, mean_vec = c(0, 1), sampCorr_mat = diag(2)),
    B = list(sampSize = 6, mean_vec = c(2, 3), sampCorr_mat = diag(2))
  )
  set.seed(2)
  out <- simulate_group_data(params, MASS::mvrnorm, "grp")

  expect_equal(nrow(out), 14)
  expect_equal(as.integer(table(out$grp)), c(8L, 6L))
})

test_that("simulate_group_data validates group definitions", {
  unnamed <- list(list(sampSize = 3, sampCorr_mat = diag(2), shape_num = c(2, 2), rate_num = c(1, 1)))
  expect_error(simulate_group_data(unnamed, generate_mvGamma_data, "grp"), "named sublists")

  missing_corr <- list(A = list(sampSize = 3, shape_num = c(2, 2), rate_num = c(1, 1)))
  expect_error(simulate_group_data(missing_corr, generate_mvGamma_data, "grp"), "sampCorr_mat")

  missing_gamma <- list(A = list(sampSize = 3, sampCorr_mat = diag(2)))
  expect_error(simulate_group_data(missing_gamma, generate_mvGamma_data, "grp"), "shape_num.*rate_num")

  missing_mean <- list(A = list(sampSize = 3, sampCorr_mat = diag(2)))
  expect_error(simulate_group_data(missing_mean, MASS::mvrnorm, "grp"), "mean_vec")

  expect_error(simulate_group_data(list(A = list(sampCorr_mat = diag(2))), stats::rnorm, "grp"), "Unsupported")
})
