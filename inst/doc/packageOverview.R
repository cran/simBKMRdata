## -----------------------------------------------------------------------------
# install.packages("simBKMRdata")


## -----------------------------------------------------------------------------
# devtools::install_github("simBKMRdata")


## -----------------------------------------------------------------------------
library(simBKMRdata)
library(tidyverse)

# Load the dataset
data("metalExposChildren_df")


analysisData_df <- metalExposChildren_df %>%
  select(QI, Cadmium:Manganese, Sex)

head(analysisData_df)

# Create a histogram with density plot for each metal
variablePlot <- analysisData_df %>%
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


## -----------------------------------------------------------------------------
param_list <- calculate_stats_gamma(
  data_df = analysisData_df, 
  group_col= "Sex", 
  using = "MoM"
)

param_list


## -----------------------------------------------------------------------------
gamma_samples <- simulate_group_data(
  param_list = param_list, 
  data_gen_fn = generate_mvGamma_data, 
  group_col_name = "Sex"
)

head(gamma_samples)


## -----------------------------------------------------------------------------
#| label: run-bkmr
library(bkmr)

set.seed(2025)

bkmrAnalysisData_df <- analysisData_df %>%
  na.omit() %>%
  mutate_at(
    vars(Cadmium:Manganese), ~ log10(. + 1) %>% as.vector
  )


# Extract the response variable (IQ) and the exposure matrix (metals)
iq <- bkmrAnalysisData_df %>% 
  pull(QI) %>% 
  na.omit()

expos <- bkmrAnalysisData_df %>%
  select(Cadmium:Manganese) %>% 
  data.matrix()

# Generate knot points using a cover design for Bayesian modeling
knots50 <- fields::cover.design(expos, nd = 50)$design

# Fit the BKMR model using MCMC
modelFit <- kmbayes(
  y = iq,         # Response variable
  Z = expos,      # Exposure matrix (metal concentrations)
  X = NULL,       # No additional covariates
  iter = 2000,    # Number of MCMC iterations
  family = "gaussian",  # Gaussian response
  verbose = TRUE,       # Output progress for each iteration
  varsel = TRUE,       # Perform variable selection
  knots = knots50      # Knot points generated earlier
)

# Extract posterior inclusion probabilities (PIPs) and sort them
pipFit <- ExtractPIPs(modelFit) %>% arrange(desc(PIP))


## -----------------------------------------------------------------------------
pipThresh_num <- calculate_pip_threshold(y = bkmrAnalysisData_df$QI)


## -----------------------------------------------------------------------------
#| echo=FALSE
variablePlot


## -----------------------------------------------------------------------------
#| echo=FALSE
library(gt)

pipFit %>%
  gt() %>%
  tab_header(
    title = "Posterior Inclusion Probabilities (PIPs)"
  )

