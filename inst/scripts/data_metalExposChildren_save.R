# Save Lucchini's Children's Metal Exposure Dataset in the Package
# Kazi Tanvir Hasan and Gabriel Odom
# 2025-04-03

# Tanvir was emailed the raw tables to create this dataset on January 15, 2024. He
# wrangled the original data set using the script <inst/scripts/data_metalExposChildren_wrangled.R>.
# He then saved the cleaned version as <"inst/extdata/metalExposChildren.csv"> on
# 2025-04-03.

# Prof. Roberto Lucchini emailed permission for us to share this data set on
# April 03, 2025. The data are shared with permission under license Creative Commons
# (CC BY 4.0), available at <https://creativecommons.org/licenses/by/4.0/>.
# Any subsequent use of the data must cite this package and acknowledge
# Prof. Roberto Lucchini (rlucchin AT fiu DOT edu).

library(tidyverse)
metalExposChildren_df <-
  read_csv("inst/extdata/metalExposChildren.csv") %>%
  select(-Zone, -Traffic) %>%
  # ADD COMMENT FROM EMAIL/DATA DICTIONARY HERE
  rename(Distance_metres = Distance, Sex = Gender) %>%
  mutate(
    SES = factor(SES, levels = c("LOW", "MEDIUM", "HIGH"), ordered = TRUE)
  ) %>%
  mutate(
    Sex = case_when(
      Sex == "F" ~ "Female", Sex == "M" ~ "Male"
    ),
    Sex = as.factor(Sex)
  )

usethis::use_data(metalExposChildren_df)
