#' @name transformers
#' @rdname transformers

#' @title Data Transformation Functions
#'
#' @description This script provides transformation functions for data scaling
#' and normalization.
#'
#' @importFrom stats sd mad
#'
NULL


#' Transform by ratio using standard deviation or MAD
#' @rdname transformers
#'
#' @param x A numeric vector or column of a dataframe to be transformed.
#' @param method Character string specifying the method: "sd" (standard
#' deviation) or "mad" (median absolute deviation).
#' @return A numeric vector where values are divided by the chosen method's
#' statistic.
#' @export
#'
#' @examples
#' trans_ratio(c(1, 2, 3, 4, 5), method = "sd")

trans_ratio <- function(x, method = c("sd", "mad")) {
  method <- match.arg(method)
  divisor <- switch(
    method,
    "sd" = sd(x, na.rm = TRUE),
    "mad" = mad(x, na.rm = TRUE)
  )
  return(x / divisor)
}


#' Transform by fractional root
#' @rdname transformers
#'
#' @param x A numeric vector or column of a dataframe to be transformed.
#' @param fracRoot The fractional power to which each element in `x` should be
#' raised. Defaults to 0.5 (square root).
#' @return A numeric vector of transformed values.
#' @export
#'
#' @examples
#' trans_root(c(1, 4, 9, 16), fracRoot = 0.5)

trans_root <- function(x, fracRoot = 0.5) {
  return(x ^ fracRoot)
}


#' Transform by logarithm with base and shift
#' @rdname transformers
#'
#' @param x A numeric vector or column of a dataframe to be transformed.
#' @param base The base of the logarithm. Defaults to 10.
#' @param shift A numeric value added to `x` before applying the logarithm to
#' avoid log(0). Defaults to 1.
#' @return A numeric vector of transformed values.
#' @export
#'
#' @examples
#' trans_log(c(1, 10, 100, 1000), base = 10, shift = 1)

trans_log <- function(x, base = 10, shift = 1) {
  return(log(x + shift, base = base))
}

