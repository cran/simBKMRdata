#' Calculate summary statistics and gamma parameters for each group
#'
#' @description This function computes the sample size, gamma distribution
#' parameters (shape and rate), and Spearman correlation matrix for each group,
#' based on the grouping column.
#'
#' @param data_df A data frame containing the data to be processed.
#' @param group_col A character string specifying the name of the column to
#' group by.
#' @param using which method will be used to estimate the multivariate Gamma 
#' shape and rate parameters. Defaults to `"MoM"` (method of moments, which was
#' used in the author's paper), or `"gMLE"` (maximum likelihood estimates from
#' the Generalized Gamma distribution without bias correction).
#'
#' @return A list of lists, where each inner list contains:
#'   - sample size (n)
#'   - sample mean vector (mean)
#'   - gamma distribution parameters (shape, rate)
#'   - Spearman correlation matrix (cor)
#' @export
#'
#' @examples
#' myData <- data.frame(
#'   GENDER = c('Male', 'Female', 'Male', 'Female', 'Male', 'Female'),
#'   VALUE1 = c(1.2, 2.3, 1.5, 2.7, 1.35, 2.5),
#'   VALUE2 = c(3.4, 4.5, 3.8, 4.2, 3.6, 4.35)
#' )
#' calculate_stats_gamma(myData, "GENDER")

calculate_stats_gamma <- function(data_df, group_col, using = c("MoM", "gMLE")) {
  
  using <- match.arg(using)
  
  # Split the data by the grouping column
  data_ls <- split(data_df, data_df[[group_col]])

  # Remove the grouping column
  data_ls <- lapply(
    X = data_ls,
    FUN = function(x) {
      x[, group_col] <- NULL
      x
    }
  )

  # Calculate statistics and moments (N, corrMat, alpha, beta)
  lapply(
    X = data_ls,
    FUN = estimate_mv_shape_rate,
    using = using
  )
  
}

