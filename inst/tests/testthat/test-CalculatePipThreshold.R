# First Unit Tests
# Kazi Tanvir Hasan and Gabriel Odom
# 2024-04-04

test_that("we can check code", {expect_true(TRUE)})

# # Test for valid input with y provided
# test_that(
#   "calculate_pip_threshold works with y input", {
#     # Test data
#     y_test <- c(2, 3, 4, 5, 6, 7)
# 
#     # Expected PIP value (this is an example, and you should verify this number)
#     expected_pip <- calculate_pip_threshold(
#       y = y_test,
#       absCV = NA,
#       sampSize = NA
#     )
# 
#     # Test if the output is numeric
#     expect_type(expected_pip, "double")
#   }
# )
# 
# # Test for valid input with absCV and sampSize
# test_that(
#   "calculate_pip_threshold works with absCV and sampSize inputs", {
#     # Test values for absCV and sample size
#     absCV_test <- 7.5
#     sampSize_test <- 300
#     
#     # Expected PIP value (this is an example, and you should verify this number)
#     expected_pip <- calculate_pip_threshold(
#       y = NA,
#       absCV = absCV_test,
#       sampSize = sampSize_test
#     )
#     
#     # Test if the output is numeric
#     expect_type(expected_pip, "double")
#   }
# )
