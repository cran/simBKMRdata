#' Helper function to estimate moment vectors/matrices for observations within
#' a group
#'
#' @param x_df A numeric data frame with observations from ONE group
#' @return A list of statistics/moments (sample size, mean, standard deviation,
#' correlation matrix, skewness) as vectors/matrices
#'
#' @importFrom stats cor
#' @export
#'
#' @examples
#' myData <- data.frame(
#'   VALUE1 = c(2.3, 2.7, 2.5),
#'   VALUE2 = c(4.5, 4.2, 4.35)
#' )
#' estimate_mv_moments(myData)
estimate_mv_moments <- function(x_df) {

  # Check if x_df is a data.frame or a matrix
  if (!is.data.frame(x_df) && !is.matrix(x_df)) {
    stop("Input x_df must be a data.frame or a matrix.")
  }

  # Convert matrix to data.frame if it is not already a data.frame
  if (is.matrix(x_df)) {
    x_df <- as.data.frame(x_df)
  }

  # Calculate sample size, mean vector, standard deviation vector, and Spearman correlation matrix
  samp_size <- nrow(x_df)
  mean_vector <- colMeans(x_df, na.rm = TRUE)  # Calculate column-wise mean, removing NA
  samp_sd <- vapply(X = x_df, FUN = sd, FUN.VALUE = numeric(1), na.rm = TRUE)

  # Spearman correlation matrix
  samp_corr <- cor(x_df, method = "spearman", use = "pairwise.complete.obs")

  # Calculate skewness for each column
  samp_skew <- vapply(
    X = seq_len(ncol(x_df)),
    FUN = function(d) {
      if (length(unique(x_df[, d])) > 1) {
        .skewness(
          x = x_df[, d],
          mean_vec = mean_vector[d],
          sampSD = samp_sd[d],
          N = samp_size
        )
      } else {
        return(NA)  # Return NA if the column has only one unique value
      }
    },
    FUN.VALUE = numeric(1)
  )

  # Return all statistics as a list
  out_ls <- list(
    sampSize = samp_size,
    mean_vec = mean_vector,
    sampSD = samp_sd,
    sampCorr_mat = samp_corr,
    sampSkew = samp_skew
  )

  return(out_ls)
}

#' Helper function to calculate skewness for a vector
#'
#' @param x A numeric vector of data
#' @param mean_vec The mean of the data
#' @param sampSD The standard deviation of the data
#' @param N The sample size
#' @return The skewness value
.skewness <- function(x, mean_vec, sampSD, N) {
  # Skewness calculation based on sample moments
  # https://en.wikipedia.org/wiki/Skewness
  numerator <- sum((x - mean_vec)^3) / N
  denominator <- (sum((x - mean_vec)^2) / (N - 1))^(3/2)
  numerator / denominator
}

