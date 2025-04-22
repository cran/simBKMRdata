#' Generate Multivariate Skewed Gamma Transformed Data
#'
#' @description This function generates multivariate normal samples, transforms
#' them into Z-scores, and then calls the `qgamma()` function to transform the
#' values for each correlated variable to those from a Gamma distribution.
#'
#' @param sampSize Number of samples to generate.
#' @param sampCorr_mat A correlation matrix for the normal distribution.
#' @param shape_num A numeric vector of shape parameters for the Gamma transformation.
#' @param rate_num A numeric vector of rate parameters for the Gamma transformation.
#' Second column: <https://en.wikipedia.org/wiki/Gamma_distribution>
#' @return A data frame containing the transformed Gamma samples.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats pnorm qgamma
#' @export
#'
#' @examples
#' p <- 4
#' N <- 1000
#' shapeGamma_num <- c(0.5, 0.75, 1, 1.25)
#' rateGamma_num <- 1:4
#' cov_mat <- diag(p)
#' generate_mvGamma_data(N, cov_mat, shapeGamma_num, rateGamma_num)

generate_mvGamma_data <- function(sampSize,
                                  sampCorr_mat,
                                  shape_num,
                                  rate_num) {

  # Check that the length of shape_num and rate_num match the number of
  # variables (p)
  p <- ncol(sampCorr_mat)
  if (length(shape_num) != p) {
    stop("The length of shape_num must match the number of variables (p).")
  }
  if (length(rate_num) != p) {
    stop("The length of rate_num must match the number of variables (p).")
  }

  # Generate multivariate normal samples
  norm_samples <- MASS::mvrnorm(
    n = sampSize,
    mu = rep.int(x = 0, times = p),
    Sigma = sampCorr_mat
  )

  # Calculate the CDF of the normal samples
  norm_cdf <- stats::pnorm(norm_samples)

  # Apply the gamma transformation
  .mvGamma_fn(norm_cdf, shape_num, rate_num)

}

# Function to apply qgamma transformation to each variable
.mvGamma_fn <- function(normCdf_mat, shape_num, rate_num) {
  dataP <- ncol(normCdf_mat)
  dataN <- nrow(normCdf_mat)

  # Apply the qgamma transformation to each column
  out_mat <- matrix(
    data = vapply(
      X = seq_len(dataP),
      FUN = function(d) {
        stats::qgamma(
          p = normCdf_mat[, d],
          shape = shape_num[d],
          rate = rate_num[d]
        )
      },
      FUN.VALUE = numeric(dataN)
    ),
    nrow = dataN, ncol = dataP, byrow = FALSE
  )

  as.data.frame(out_mat)

}
