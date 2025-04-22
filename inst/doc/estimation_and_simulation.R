## ----setup, include=TRUE------------------------------------------------------
# Load necessary libraries
library(MASS)         # For generating multivariate normal data
library(simBKMRdata)  # For generating skewed Gamma data and estimating moments


## -----------------------------------------------------------------------------
# Example using MASS::mvrnorm for normal distribution
param_list <- list(
  Group1 = list(
    mean_vec = c(1, 2), 
    sampCorr_mat = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
    sampSize = 100
  ),
  Group2 = list(
    mean_vec = c(2, 3), 
    sampCorr_mat = matrix(c(1, 0.3, 0.3, 1), 2, 2), 
    sampSize = 150
  )
)

mvnorm_samples <- simulate_group_data(param_list, MASS::mvrnorm, "Group")


## -----------------------------------------------------------------------------
# Plot the first two variables of the multivariate normal data
plot(
  mvnorm_samples[, 1], mvnorm_samples[, 2],
  main = "Scatterplot: MV Normal Data",
  xlab = "Variable 1", ylab = "Variable 2",
  pch = 19, col = "blue"
)


## -----------------------------------------------------------------------------
myData <- data.frame(
  GENDER = c('Male', 'Female', 'Male', 'Female', 'Male', 'Female'),
  VALUE1 = c(1.2, 2.3, 1.5, 2.7, 1.35, 2.5),
  VALUE2 = c(3.4, 4.5, 3.8, 4.2, 3.6, 4.35)
)
calculate_stats_gaussian(data_df = myData, group_col = "GENDER")


## -----------------------------------------------------------------------------
# Example using generate_mvGamma_data for Gamma distribution
param_list <- list(
   Group1 = list(
     sampCorr_mat = matrix(c(1, 0.5, 0.5, 1), 2, 2),
     shape_num = c(2, 2), 
     rate_num = c(1, 1), 
     sampSize = 100
   ),
   Group2 = list(
     sampCorr_mat = matrix(c(1, 0.3, 0.3, 1), 2, 2),
     shape_num = c(2, 2), 
     rate_num = c(1, 1), 
     sampSize = 150
   )
 
)

gamma_samples <- simulate_group_data(
  param_list, generate_mvGamma_data, "Group"
)


## -----------------------------------------------------------------------------
# Plot the density of the first and second variable for Gamma data
old_par_mfrow <- par()[["mfrow"]]
par(mfrow = c(2, 1))
plot(density(gamma_samples[, 1]), main = "Gamma Variable 1", col = "blue")
plot(density(gamma_samples[, 2]), main = "Gamma Variable 2", col = "blue")
par(mfrow = old_par_mfrow)


## -----------------------------------------------------------------------------
myData <- data.frame(
   GENDER = c('Male', 'Female', 'Male', 'Female', 'Male', 'Female'),
   VALUE1 = c(1.2, 2.3, 1.5, 2.7, 1.35, 2.5),
   VALUE2 = c(3.4, 4.5, 3.8, 4.2, 3.6, 4.35)
)
calculate_stats_gamma(data_df = myData, group_col= "GENDER", using = "MoM")

