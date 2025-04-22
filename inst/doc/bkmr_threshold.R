## -----------------------------------------------------------------------------
# Load the necessary libraries
library(bkmr)
library(simBKMRdata)
library(tidyverse)


## -----------------------------------------------------------------------------
# Set the seed for reproducibility
set.seed(2025)

# Load the dataset and preprocess
data("metalExposChildren_df")

bkmrAnalysisData_df <- 
  metalExposChildren_df %>%
  select(QI, Cadmium:Manganese) %>%  # Selecting relevant columns
  na.omit() %>%  # Remove rows with missing values
  mutate_at(
    vars(Cadmium:Manganese), ~ log10(. + 1) %>% as.vector
  )  # Log-transform metal concentrations


## -----------------------------------------------------------------------------
# Create a histogram with density plot for each metal
metalDataPlot <- bkmrAnalysisData_df %>%
  select(Cadmium:Manganese) %>%
  pivot_longer(cols = everything(), names_to = "Metal", values_to = "Value") %>%
  ggplot(aes(x = Value, fill = Metal)) + 
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.5, position = "identity") + 
  geom_density(aes(color = Metal), linewidth = 1) + 
  facet_wrap(~Metal, scales = "free") + 
  theme_minimal() + 
  labs(title = "Histogram with Density Plot for Metals", 
       x = "Concentration", 
       y = "Density") + 
  theme(legend.position = "none")

# Display the plot
metalDataPlot


## -----------------------------------------------------------------------------
# Extract the response variable (IQ score or QI)
iq <- bkmrAnalysisData_df %>%
  pull(QI) %>%
  na.omit()

# Convert exposure variables (metals) to a matrix for modeling
expos <- data.matrix(bkmrAnalysisData_df %>% select(Cadmium:Manganese))


## -----------------------------------------------------------------------------
# Generate knot points using a cover design for Bayesian modeling
knots50 <- fields::cover.design(expos, nd = 50)$design


## -----------------------------------------------------------------------------
#| label: run-BKMR-MCMC
# Fit the BKMR model using MCMC
modelFit <- kmbayes(
  y = iq,         # Response variable
  Z = expos,      # Exposure matrix (metal concentrations)
  X = NULL,       # No additional covariates
  iter = 1000,    # Number of MCMC iterations
  family = "gaussian",  # Gaussian response
  verbose = TRUE,       # Output progress for each iteration
  varsel = TRUE,       # Perform variable selection
  knots = knots50      # Knot points generated earlier
)


## -----------------------------------------------------------------------------
# Extract posterior inclusion probabilities (PIPs) and sort them
pipFit <- ExtractPIPs(modelFit) %>%
  arrange(desc(PIP))

# Display the PIPs
pipFit


## -----------------------------------------------------------------------------
# Calculate the dynamic threshold for model inclusion
pipThresh_fn <- calculate_pip_threshold(
  absCV = sd(bkmrAnalysisData_df$QI) / mean(bkmrAnalysisData_df$QI),  # Coefficient of variation
  sampSize = length(bkmrAnalysisData_df$QI)  # Sample size
)

# Display the threshold
pipThresh_fn


## -----------------------------------------------------------------------------
# Identify exposures with PIPs greater than the threshold
significant_exposures <- pipFit %>%
  filter(PIP > pipThresh_fn)

# Display significant exposures
significant_exposures

