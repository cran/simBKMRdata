#' Simulate Group Data
#'
#' @description This function generates data for each group by invoking the
#' specified data generation function once per group. It binds the generated
#' data together into a single data frame.
#'
#' @param param_list A list of named sublists, where each sublist contains the
#' parameters for a group (mean, shape, rate, etc.). The list must be named with
#' group names that match the groupings stated in `group_col_name`.
#' @param data_gen_fn A function for data generation. Currently we can choose
#' either `generate_mvGamma_data` or `MASS::mvrnorm`.
#' @param group_col_name The name of the grouping/label column to be created in
#' the final data frame.
#'
#' @return A data frame with the simulated data for all groups, including the grouping column.
#'
#' @examples
#' # Example using MASS::mvrnorm for normal distribution
#' param_list <- list(
#'   Group1 = list(mean_vec = c(1, 2), sampCorr_mat = matrix(c(1, 0.5, 0.5, 1), 2, 2), sampSize = 100),
#'   Group2 = list(mean_vec = c(2, 3), sampCorr_mat = matrix(c(1, 0.3, 0.3, 1), 2, 2), sampSize = 150)
#' )
#' simulate_group_data(param_list, MASS::mvrnorm, "Group")
#'
#' # Example using generate_mvGamma_data for Gamma distribution
#' param_list <- list(
#'   Group1 = list(sampCorr_mat = matrix(c(1, 0.5, 0.5, 1), 2, 2),
#'                 shape_num = c(2, 2), rate_num = c(1, 1), sampSize = 100),
#'   Group2 = list(sampCorr_mat = matrix(c(1, 0.3, 0.3, 1), 2, 2),
#'                 shape_num = c(2, 2), rate_num = c(1, 1), sampSize = 150)
#' )
#' simulate_group_data(param_list, generate_mvGamma_data, "Group")
#'
#' @export
simulate_group_data <- function(param_list, data_gen_fn, group_col_name) {

  # Check if the list of parameters is named
  if (is.null(names(param_list))) {
    stop("param_list must be a list of named sublists for each group.")
  }

  # Initialize an empty list to store the generated data for each group
  all_data <- list()

  # Iterate over each group in param_list
  for (group_name in names(param_list)) {

    # Get the parameters for the current group
    group_params <- param_list[[group_name]]

    # Check if the required parameters are present in the group
    if (is.null(group_params$sampCorr_mat)) {
      stop("Each group must have 'sampCorr_mat' in its parameters.")
    }

    # Check for distribution-specific parameters and generate data
    if (identical(data_gen_fn, generate_mvGamma_data)) {
      # Ensure Gamma-specific parameters are present
      if (is.null(group_params$shape_num) || is.null(group_params$rate_num)) {
        stop("Each group must have 'shape_num' and 'rate_num' for the Gamma distribution.")
      }

      # Generate data using Gamma distribution
      generated_data <- data_gen_fn(
        sampSize = group_params$sampSize,
        sampCorr_mat = group_params$sampCorr_mat,
        shape_num = group_params$shape_num,
        rate_num = group_params$rate_num
      )

    } else if (identical(data_gen_fn, MASS::mvrnorm)) {
      # Generate data using Normal distribution
      if (is.null(group_params$mean_vec)) {
        stop("Each group must have 'mean_vec' in its parameters.")
      }

      generated_data <- data_gen_fn(
        n = group_params$sampSize,
        mu = group_params$mean_vec,
        Sigma = group_params$sampCorr_mat
      )

    } else {
      stop("Unsupported data generation function provided.")
    }

    # Create a data frame for the generated data and add the group label
    generated_data_df <- as.data.frame(generated_data)
    generated_data_df[[group_col_name]] <- group_name  # Add group label

    # Append the data to the list
    all_data[[group_name]] <- generated_data_df
  }

  # Combine all the group data frames into one
  combined_data <- do.call(rbind, all_data)
  # currently it's adding non-unique row name = group label
  # (female.1, female.2,...)
  rownames(combined_data) <- NULL

  return(combined_data)
}
