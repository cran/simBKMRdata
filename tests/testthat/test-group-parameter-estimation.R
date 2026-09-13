test_that("calculate_stats_gaussian estimates each group separately", {
  x <- data.frame(
    grp = rep(c("A", "B"), each = 4),
    x1 = c(1, 2, 3, 4, 10, 11, 12, 13),
    x2 = c(2, 3, 5, 8, 4, 6, 9, 13)
  )
  out <- calculate_stats_gaussian(x, "grp")

  expect_named(out, c("A", "B"))
  expect_equal(out$A$sampSize, 4)
  expect_equal(out$B$sampSize, 4)
  expect_equal(unname(out$A$mean_vec), c(2.5, 4.5))
})

test_that("calculate_stats_gaussian validates the group column", {
  x <- data.frame(group = c("A", "B"), x = c(1, 2))
  expect_error(calculate_stats_gaussian(x, "missing"), "Grouping column not found")
})

test_that("calculate_stats_gamma estimates each group separately", {
  x <- data.frame(
    grp = rep(c("A", "B"), each = 5),
    x1 = c(1, 2, 3, 4, 5, 2, 3, 4, 5, 6),
    x2 = c(2, 3, 5, 7, 11, 1, 2, 4, 8, 16)
  )
  out <- calculate_stats_gamma(x, "grp", using = "MoM")

  expect_named(out, c("A", "B"))
  expect_equal(out$A$sampSize, 5)
  expect_length(out$A$shape_num, 2)
  expect_true(all(out$A$shape_num > 0))
  expect_true(all(out$A$rate_num > 0))
})

test_that("calculate_stats_gamma validates estimation method", {
  x <- data.frame(grp = rep(c("A", "B"), each = 2), x = 1:4)
  expect_error(calculate_stats_gamma(x, "grp", using = "bad"))
})
