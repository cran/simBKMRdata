# Wrangled Lucchini's Children's Metal Exposure Dataset in the Package
# Kazi Tanvir Hasan and Gabriel Odom
# 2025-04-03


# # Tanvir (Me) was emailed the raw tables to create this dataset on
# January 15, 2024. In this script I wrangled the original data set using this
# script saved the cleaned version as "inst/extdata/metalExposChildren.csv"

library(readxl)
library(tidyverse)

set.seed(2025) # Ensure reproducibility

##### CCM-ISEIA common/overlapping data ######

# WISC-IV

wiscIV_CCM <- read_excel("data/WISC_TA.xlsx") %>%
  select(ID, age, QI)

wiscIV_ISEA <- read_excel("data/WISC IV_ISEA_sharing.xlsx", sheet = "data") %>%
  select(ID, age, QI) %>%
  # Exclude data collected in Brescia
  filter(ID %in% paste0("TA", 314:632))  # Creates a character vector from TA314 to TA632

wiscIV_df <- rbind(wiscIV_CCM, wiscIV_ISEA)


# Standard Progressive Matrices (SPM) parent

spm_CCM <- read_excel("data/SRS-SPM_TA.xlsx") %>%
  select(
    ID, SPMT_cor
  ) %>%
  rename(SPM = SPMT_cor)

spm_ISEA <- read_delim(
  "docs/ISEIA dataset/ISEIAALL_DATA_2023-10-12_1435_sharing.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
) %>%
  select(
    raven_subject_id, spm_totesatte
  ) %>%
  rename(
    ID = raven_subject_id,
    SPM = spm_totesatte
  ) %>%
  filter(ID %in% paste0("TA", 314:632))  # Creates a character vector from TA314 to TA632

spm_df <- rbind(spm_CCM, spm_ISEA)

# Home questionnaire

socioDemo_CCM <- read_excel("data/SOCIODEMO_TA.xlsx") %>%
  select(
    ID, Gender, BMI, Children_SES_index, Zone, Amount_of_traffic, Distance
  ) %>%
  rename(
    SES = Children_SES_index,
    Traffic = Amount_of_traffic
  )


socioDemo1_ISEA <- read_delim(
  "data/QI_soil_dataset.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE
) %>%
  select(ID, sex, SES, zone, distance) %>%
  rename(
    Gender = sex,
    Zone = zone,
    Distance = distance
  )

socioDemo2_ISEA <- ISEIA_dataset <- read_excel("data/ISEIA_dataset.xlsx") %>%
  select(cbcl_teacher_subject_id,weight, height, via74_2) %>%
  rename(
    ID = cbcl_teacher_subject_id,
    Traffic = via74_2
  )

socioDemo_ISEA <- left_join(socioDemo1_ISEA, socioDemo2_ISEA, by = "ID") %>%
  mutate(
    BMI = weight/(height * 0.01)^2
  ) %>%
  select(ID, Gender, BMI, SES, Zone, Traffic, Distance)


socioDemo_df <- rbind(socioDemo_CCM, socioDemo_ISEA)

# Blood metals

melats_CCM <- read_excel("data/BIOMARKERS_TA.xlsx", sheet = "bio2") %>%
  rename(
    Cadmium = "Cd-U",
    Mercury = "Hg-U",
    Arsenic = "As-U",
    Lead = "Pb-B",
    Manganese = "Mn-H"
  ) %>%
  select(ID, Cadmium, Mercury, Arsenic, Lead, Manganese)


metals_ISEA <- read_excel("data/Blood_HeavyMetals_sharing.xlsx") %>%
  rename(ID= "Sample ID")

manganese_ISEA <- read_delim(
  "data/QI_soil_dataset.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE
) %>%
  select(ID, Mn) %>%
  rename(Manganese = Mn)

metalsComplete_ISEA <- left_join(metals_ISEA, manganese_ISEA, by = "ID") %>%
  select(ID, Cadmium, Mercury, Arsenic, Lead, Manganese)

metals_df <- rbind(melats_CCM, metalsComplete_ISEA)
# Teeth
teethMetals_df <- read_excel("data/RL_TEETH_20250122_tp_bpj.xlsx") %>%
  slice(-1) %>%
  select(id, sample, spot, nl, Cd, Hg, As, Pb, Mn) %>%
  rename(
    ID = id,
    Cadmium = Cd,
    Mercury = Hg,
    Arsenic = As,
    Lead = Pb,
    Manganese = Mn
  ) %>%
  mutate(
    across(
      c(Cadmium, Mercury, Arsenic, Lead, Manganese),
      ~ as.numeric(.)
    )
  ) %>%
  mutate(ID_tooth = str_extract(string = ID, pattern = "A$|B$")) %>%
  mutate(ID = sub("A$|B$", "", ID))

# Analysis Data

analysisData_df <- wiscIV_df %>%
  left_join(., metals_df, by = "ID") %>%
  left_join(., socioDemo_df, by = "ID" ) %>%
  left_join(., spm_df, by = "ID")

analysisData_df <- analysisData_df %>%
  mutate(
    Gender = case_when(
      Gender == "0" ~ "M",
      Gender == "1" ~ "F",
      TRUE ~ Gender
    ) %>% as.factor()
  ) %>%
  mutate(
    SES = case_when(
      SES == "0" ~ "LOW",
      SES == "1" ~ "MEDIUM",
      SES == "2" ~ "HIGH",
      TRUE ~ SES
    ) %>% as.factor()
  ) %>%
  mutate(
    Zone = as.factor(Zone),
    Traffic = as.factor(Traffic)
  ) %>%
  mutate(
    across(
      -c(ID, Gender, SES, Zone, Traffic),
      ~ as.numeric(as.character(.))
    )  # Convert all remaining columns to numeric
  )

# Data Partition

bkmrData_df <- analysisData_df %>%
  filter(
    !(ID %in% unique(teethMetals_df$ID))
  )


# # Save bkmrData_df as an csv file
# write_csv(bkmrData_df, path = "inst/extdata/metalExposChildren.csv")

