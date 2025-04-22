#' Simulate Group Multivariate Gaussian Data
#'
#' @description This function generates data for each group from a Multivariate
#' Gaussian (Normal) Distribution by invoking this distribution's random
#' generator once per group. It binds the generated data together into a single
#' data frame.
#'
#' @param param_list A list of named sublists, where each sublist contains the
#' parameters for a group (sample size, mean, standard deviation, and
#' correlation matrix). The dimension of the parameters for each group must be
#' the same.
#' @param group_col_name The column name of the grouping/label column to be
#' created in the final data frame. The values are taken from the names of the
#' sublists of `param_list`. Defaults to "group". See the example below.
#'
#' @return A data frame with the simulated data for all groups, including the
#' grouping column.
#' 
#' 
#' @export
#' @importFrom MASS mvrnorm
#'
#' @examples
#' # Example using MASS::mvrnorm for normal distribution
#' param_list <- list(
#'   Male = list(
#'     sampSize = 50,
#'     mean_vec = c(1, 2),
#'     sampSD = c(2, 1),
#'     sampCorr_mat = matrix(c(1, 0.5, 0.5, 1), 2, 2)
#'   ),
#'   Female = list(
#'     sampSize = 100,
#'     mean_vec = c(2, 3),
#'     sampSD = c(1, 2),
#'     sampCorr_mat = matrix(c(1, 0.3, 0.3, 1), 2, 2)
#'   )
#' )
#' simulate_group_gaussian(param_list, "Sex")
#' 
simulate_group_gaussian <- function(param_list, group_col_name) {

  # Check if the list of parameters is named
  group_names_all <- names(param_list)
  if (is.null(group_names_all)) {
    stop(
      "param_list must be a list of named sublists for each group. The names
      will be used as the values for the group_col_name column.")
  }
  
  # Add check that all the dimensions are equal across the groups
  
  # Expected parameter names
  gaussian_params_expected <- c(
    "sampSize", "mean_vec", "sampSD", "sampCorr_mat"
  )
  
  simOut_ls <- lapply(
    X = group_names_all,
    FUN = function(group_name) {
      
      # Get the parameters for the current group
      group_params <- param_list[[group_name]]
      param_names <- names(group_params)
      param_check_lgl <- gaussian_params_expected %in% param_names
      if (!all(param_check_lgl)) {
        stop(
          "All parameter sublists must contain sample size (sampSize), sample mean
        (mean_vec), sample standard deviation (sampSD), and sample correlation
        matrix (sampCorr_mat); these list names must be EXACT.")
      }
      
      generated_data <- MASS::mvrnorm(
        n = group_params$sampSize,
        mu = group_params$mean_vec,
        Sigma = group_params$sampCorr_mat
      )
      
      # Create a data frame for the generated data and add the group label
      generated_data_df <- as.data.frame(generated_data)
      generated_data_df[[group_col_name]] <- group_name  # Add group label
      
      generated_data_df
      
    }
  )

  # Combine all the group data frames into one
  do.call("rbind", simOut_ls)
  
}

# # Initialize an empty list to store the generated data for each group
# all_data <- list()
# 
# # Iterate over each group in param_list
# for (group_name in names(param_list)) {
#   
#   # Get the parameters for the current group
#   group_params <- param_list[[group_name]]
#   gaussian_param_names <- names(group_params)
#   gaussian_params_expected <- c(
#     "sampSize", "mean_vec", "sampSD", "sampCorr_mat"
#   )
#   param_check_lgl <- gaussian_params_expected %in% gaussian_param_names
#   if (!all(param_check_lgl)) {
#     stop(
#       "All parameter sublists must contain sample size (sampSize), sample mean
#         (mean_vec), sample standard deviation (sampSD), and sample correlation
#         matrix (sampCorr_mat); these list names must be EXACT.")
#   }
#   
#   # # Check if the required parameters are present in the group
#   # if (is.null(group_params$mean_vec) || is.null(group_params$sampCorr_mat)) {
#   #   stop("Each group must have 'mean_vec' and 'sampCorr_mat' in its parameters.")
#   # }
#   
#   generated_data <- MASS::mvrnorm(
#     n = group_params$sampSize,
#     mu = group_params$mean_vec,
#     Sigma = group_params$sampCorr_mat
#   )
#   
#   # Create a data frame for the generated data and add the group label
#   generated_data_df <- as.data.frame(generated_data)
#   generated_data_df[[group_col_name]] <- group_name  # Add group label
#   
#   # Append the data to the list
#   all_data[[group_name]] <- generated_data_df
# }

