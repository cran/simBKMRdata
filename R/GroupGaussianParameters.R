#' Calculate summary statistics for each group
#'
#' @description This function computes the sample size, mean vector, standard
#' deviation vector, Spearman correlation matrix, and skewness vector for each
#' group, based on the grouping column.
#'
#' @param data_df A data frame containing the data to be processed.
#' @param group_col A character string specifying the name of the column to
#' group by.
#'
#' @return A list of lists, where each inner list contains the following
#' parameter estimates for one group:
#'   - sample size (`sampSize`)
#'   - sample mean vector (`xBar`)
#'   - sample standard deviation vector (`sampSD`)
#'   - sample Spearman correlation matrix (`sampCorr_mat`)
#'   - sample skewness (`sampSkew`)
#' @export
#'
#' @examples
#' myData <- data.frame(
#'   GENDER = c('Male', 'Female', 'Male', 'Female', 'Male', 'Female'),
#'   VALUE1 = c(1.2, 2.3, 1.5, 2.7, 1.35, 2.5),
#'   VALUE2 = c(3.4, 4.5, 3.8, 4.2, 3.6, 4.35)
#' )
#' calculate_stats_gaussian(myData, "GENDER")

calculate_stats_gaussian <- function(data_df, group_col) {
  # Check if the grouping column exists
  if (!group_col %in% names(data_df)) {
    stop("Grouping column not found in the data frame.")
  }

  # Split the data by the grouping column
  data_ls <- split(data_df, data_df[[group_col]])

  # Remove the grouping column and apply statistics calculation for each group
  data_ls <- lapply(
    X = data_ls,
    FUN = function(x) {
      x[, group_col] <- NULL
      x
    }
  )

  # Calculate statistics and moments (N, xBar, sd, corrMat, skew)
  lapply(
    X = data_ls,
    FUN = estimate_mv_moments
  )
  
}
