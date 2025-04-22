# check all the functions. if there is matrix change it to data.frame
#
# we are estimating 1/beta

param_list <- list(
  Group1 = list(sampCorr_mat = matrix(c(1, 0, 0, 1), 2, 2),
                shape_num = c(10, 6), rate_num = c(5, 2), sampSize = 10000),
  Group2 = list(sampCorr_mat = matrix(c(1, 0.3, 0.3, 1), 2, 2),
                shape_num = c(2, 2), rate_num = c(1, 1), sampSize = 10)
)

example_df <- simulate_group_data(param_list, generate_mvGamma_data, "Group")

example_df <- example_df %>%
  filter(
    Group == "Group1"
  )

example_df %>%
  summarise(
    # should = 10/5 = 2 if rate; 10 * 5 = 50 if scale
    meanV1 = mean(V1),
    # should = 10/5^2 = 0.4 if rate; 10 * 5^2 = 250 if scale
    varV1 = var(V1),
    # should = 6/2 = 3 if rate; 6 * 2 = 12 if scale
    meanV2 = mean(V2),
    # should = 6/2^2 = 1.5 if rate; 6 * 2^2 = 24 if scale
    varV2 = var(V2),
  )



example_df %>%
  select(-Group) %>%
  estimate_mv_shape_rate(using = "MoM")
# should be identity correlation matrix, shape vector close to (10, 6),
# rate vector close to (5, 2)

example_df %>%
  select(-Group) %>%
  estimate_mv_shape_rate(using = "gMLE")
# should be identity correlation matrix, shape vector close to (10, 6),
# rate vector close to (5, 2)

