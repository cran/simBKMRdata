## -----------------------------------------------------------------------------
#| eval: false
# install.packages("simBKMRdata")


## -----------------------------------------------------------------------------
#| eval: false
# devtools::install_github("khasa006/simBKMRdata")


## -----------------------------------------------------------------------------
#| message: false
library(tidyverse)
library(simBKMRdata)
library(gt)


## -----------------------------------------------------------------------------
data("metalExposChildren_df")

analysisData_df <- metalExposChildren_df %>%
  select(QI, Cadmium:Manganese, Sex)

head(analysisData_df) %>%
  gt() %>%
  tab_header(
    title = "Top Rows of the Data"
  ) %>%
  fmt_number(
    decimals = 2
  ) %>%
  opt_table_outline()


## -----------------------------------------------------------------------------
#| label: fig-metals-density
#| fig-cap: Histogram with Density Plot for Metals
#| warning: false

ggplot(
  data = analysisData_df %>%
    select(Cadmium:Manganese) %>%
    pivot_longer(
      cols = everything(), names_to = "Metal", values_to = "Value"
    )
) + 
  aes(x = Value, fill = Metal) + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  labs(
    title = "Histogram with Density Plot for Metals", 
    x = "Concentration", 
    y = "Density"
  ) +
  geom_histogram(
    aes(y = ..density..), bins = 30, alpha = 0.5, position = "identity"
  ) + 
  geom_density(aes(color = Metal), linewidth = 1) + 
  facet_wrap(~Metal, scales = "free")


## -----------------------------------------------------------------------------
param_list <- calculate_stats_gamma(
  data_df = analysisData_df, 
  group_col = "Sex", 
  using = "MoM"
)

# Examine list contents
str(
  param_list,
  max.level = 2,
  give.attr = FALSE
)


## -----------------------------------------------------------------------------
gamma_samples <- simulate_group_data(
  param_list = param_list, 
  data_gen_fn = generate_mvGamma_data, 
  group_col_name = "Sex"
)

head(gamma_samples) %>% 
  gt() %>%
  tab_header(
    title = "Top Rows of the gamma Sample"
  ) %>%
  fmt_number(
    decimals = 2
  ) %>%
  opt_table_outline()


## -----------------------------------------------------------------------------
pipThresh_num <- calculate_pip_threshold(
  y = na.omit(analysisData_df$QI)
)
pipThresh_num


## -----------------------------------------------------------------------------
#| label: transform-data-for-bkmr

bkmrAnalysisData_df <- 
  analysisData_df %>%
  drop_na() %>%
  mutate_at(
    vars(Cadmium:Manganese), ~ trans_log(.) %>% as.vector
  )


## -----------------------------------------------------------------------------
#| label: create-bkmr-data-objects
set.seed(2025)
expos <- bkmrAnalysisData_df %>%
  select(Cadmium:Manganese) %>% 
  data.matrix()
is_female <- (bkmrAnalysisData_df$Sex == "Female") + 0L

# Generate knot points using a cover design for Bayesian modeling
knots50 <- fields::cover.design(expos, nd = 50)$design


## -----------------------------------------------------------------------------
#| label: run-bkmr
#| message: false
library(bkmr)

set.seed(2025)

# Fit the BKMR model using MCMC
modelFit <- kmbayes(
  # Response variable
  y = bkmrAnalysisData_df$QI,
  # Exposure matrix (metal concentrations)
  Z = expos,     
  # Sex at birth covariate
  X = is_female,   
  # Number of MCMC iterations; set to at least 10000 for real analysis
  iter = 2000,         
  family = "gaussian", # Gaussian response
  verbose = FALSE,     # Output progress for each iteration
  varsel = TRUE,       # Perform variable selection
  knots = knots50      # Knot points generated earlier
)


## -----------------------------------------------------------------------------
#| label: fig-Univariate-exposure-response-function
#| warning: false
# Generate univariate predictor-response relationships
predRespUnivar <- PredictorResponseUnivar(fit = modelFit)

# Plot univariate predictor-response functions
predRespUnivarPlot <- ggplot(
  predRespUnivar, 
  aes(
    z, 
    est, 
    ymin = est - 1.96 * se, 
    ymax = est + 1.96 * se
  )
) + 
  geom_smooth(stat = "identity") +  # Add smooth lines with confidence intervals
  facet_wrap(~ variable) +          # Create separate plots for each variable
  ylab("h(z)")                      # Label the y-axis



## -----------------------------------------------------------------------------
#| label: tbl-PIPs
#| tbl-cap: Posterior Inclusion Probabilities (PIPs) for Metal Exposures
#| echo: false

pipFit <- ExtractPIPs(modelFit) %>% arrange(desc(PIP))
pipFit %>%
  gt() %>%
  tab_header(
    title = "Posterior Inclusion Probabilities (PIPs)"
  )


## -----------------------------------------------------------------------------
#| label: fig-Univariate-exposure-response-plot
#| fig-cap:  Univariate exposure-response function plot
#| warning: false

# Plot univariate predictor-response functions
predRespUnivarPlot <- ggplot(
  predRespUnivar, 
  aes(
    z, 
    est, 
    ymin = est - 1.96 * se, 
    ymax = est + 1.96 * se
  )
) + 
  geom_smooth(stat = "identity") +  # Add smooth lines with confidence intervals
  facet_wrap(~ variable) +          # Create separate plots for each variable
  ylab("h(z)")                      # Label the y-axis

predRespUnivarPlot

