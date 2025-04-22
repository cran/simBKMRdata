#' Simulate Group Multivariate Data
#'
#' @description This function generates data for each group from a Multivariate
#' Gamma Distribution by invoking this distribution's random generator once
#' per group. It binds the generated data together into a single data frame.
#'
#' @param param_list A list of named sublists, where each sublist contains the
#' parameters for a group (sample size, mean, standard correlation matrix,
#' shape, and rate parameter). The dimension of the parameters for each group
#' must be the same.
#' @param group_col_name The column name of the grouping/label column to be
#' created in the final data frame. The values are taken from the names of the
#' sublists of `param_list`. Defaults to "group". See the example below.
#'
#' @return A data frame with the simulated data for all groups, including the
#' grouping column.
#'
#' @export
#'
#' @examples
#'
#' # Example using generate_mvGamma_data for MV Gamma distribution
#' param_list <- list(
#'  Male = list(
#'   sampSize = 100,
#'   sampCorr_mat = matrix(c(1, 0.5, 0.5, 1), 2, 2),  # Covariance matrix
#'   shape_num = c(2, 2),  # Shape parameters for Gamma distribution
#'   rate_num = c(1, 1)    # Rate parameters for Gamma distribution
#'  ),
#'  Female = list(
#'   sampSize = 150,
#'   sampCorr_mat = matrix(c(1, 0.3, 0.3, 1), 2, 2),
#'   shape_num = c(1, 4),
#'   rate_num = c(0.5, 2)
#'  )
#' )
#' simulate_group_gamma(param_list, "Sex")
#'

simulate_group_gamma <- function(param_list, group_col_name) {

  # Check if the list of parameters is named
  group_names_all <- names(param_list)
  if (is.null(group_names_all)) {
    stop(
      "param_list must be a list of named sublists for each group. The names
      will be used as the values for the group_col_name column.")
  }

  # Add check that all the dimensions are equal across the groups

  # Expected parameter names
  gamma_params_expected <- c(
    "sampSize", "sampCorr_mat", "shape_num", "rate_num"
  )

  simOut_ls <- lapply(
    X = group_names_all,
    FUN = function(group_name) {

      # Get the parameters for the current group
      group_params <- param_list[[group_name]]
      param_names <- names(group_params)
      param_check_lgl <- gamma_params_expected %in% param_names
      if (!all(param_check_lgl)) {
        stop(

          "All parameter sublists must contain sample size (sampSize), sample
          correlation matrix (sampCorr_mat), sample shape parameter (shape_num)
          and sample rate parameter (rate_num); these list names must be EXACT."
        )
      }

      generated_data <- generate_mvGamma_data(
        sampSize = group_params$sampSize,  # Number of samples
        sampCorr_mat = group_params$sampCorr_mat,  # Pre-estimated cov matrix
        shape_num = group_params$shape_num,  # Pre-estimated shape parameters
        rate_num = group_params$rate_num  # Pre-estimated rate parameters
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




# simulate_group_gamma <- function(param_list, group_col_name) {
#
#   # Check if the list of parameters is named
#   if (is.null(names(param_list))) {
#     stop("param_list must be a list of named sublists for each group.")
#   }
#
#   # Initialize an empty list to store the generated data for each group
#   all_data <- list()
#
#   # Iterate over each group in param_list
#   for (group_name in names(param_list)) {
#
#     # Get the parameters for the current group
#     group_params <- param_list[[group_name]]
#
#     # Check if the required parameters are present in the group
#     if (
#       is.null(group_params$mean_vec) || is.null(group_params$sampCorr_mat)
#     ) {
#       stop(
#         "Each group must have 'mean_vec' and 'sampCorr_mat' in its
#         parameters."
#       )
#     }
#
#     # Check for distribution-specific parameters and generate data
#     if (identical(data_gen_fn, generate_mvGamma_data)) {
#       # Ensure Gamma-specific parameters are present
#       if (is.null(group_params$shape_num) || is.null(group_params$rate_num)) {
#         stop(
#           "Each group must have 'shape_num' and 'rate_num' for the Gamma
#           distribution."
#         )
#       }
#
#       # Generate data using Gamma distribution
#       generated_data <- data_gen_fn(
#         sampSize = group_params$sampSize,
#         mean_vec = group_params$mean_vec,
#         sampCorr_mat = group_params$sampCorr_mat,
#         shape_num = group_params$shape_num,
#         rate_num = group_params$rate_num
#       )
#
#     } else if (identical(data_gen_fn, MASS::mvrnorm)) {
#       # Generate data using Normal distribution
#       generated_data <- data_gen_fn(
#         n = group_params$sampSize,
#         mu = group_params$mean_vec,
#         Sigma = group_params$sampCorr_mat
#       )
#
#     } else {
#       stop("Unsupported data generation function provided.")
#     }
#
#     # Create a data frame for the generated data and add the group label
#     generated_data_df <- as.data.frame(generated_data)
#     generated_data_df[[group_col_name]] <- group_name  # Add group label
#
#     # Append the data to the list
#     all_data[[group_name]] <- generated_data_df
#   }
#
#   # Combine all the group data frames into one
#   combined_data <- do.call(rbind, all_data)
#
#   return(combined_data)
# }
