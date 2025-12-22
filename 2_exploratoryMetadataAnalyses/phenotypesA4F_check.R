library(tidyverse)
library(janitor)
library(readxl)
library(reshape2)

rm(list = ls())

# dataDir is a clone of AI4Food's GitHub
dataDir <- paste0(iniDir, "AI4FoodDB/")

# DATA IMPORT ##################################################################
setwd(dataDir)

participant <- read_csv("./participant_information.csv") 

vol_ids <- participant %>% pull(id)
men_ids <- participant %>% filter(sex == "Male") %>% pull(id)
women_ids <- participant %>% filter(sex == "Female") %>% pull(id)


anthropo_data <- read_csv('./DS1_AnthropometricMeasurements/anthropometrics.csv') %>%
  filter(visit != 2) %>%
  left_join(participant %>% select(id, age), 
            by = "id")

vitalsigns_data <- read_csv("DS6_VitalSigns/vital_signs.csv")

bq_data <- read_csv('./DS4_Biomarkers/biomarkers.csv') %>%
  mutate(crp_mg_dl = if_else(crp_mg_dl == "<0.003", NA, 
                             as.numeric(crp_mg_dl)),
         hba1c_perc = 100 * hba1c_perc)

med_data <- read_csv("./DS2_LifestyleHealth/health.csv")


# SELECT VARIABLES OF INTEREST ------
d_vol <- participant %>%
  select(id, sex, age)

d_anthropo <- anthropo_data %>%
  select(id, visit, bmi_kg_m2, waist_cm, hip_cm) %>%
  filter(visit != 2) %>%
  mutate(visit = if_else(visit == 0, 1, visit))

d_vs <- vitalsigns_data %>%
  filter(visit != 2) %>%
  mutate(visit = if_else(visit == 0, 1, visit))

d_bq <- bq_data %>%
  select(id, visit,
         hba1c_perc, insulin_uui_ml, homa, glu_mg_dl, 
         chol_mg_dl, tri_mg_dl, hdl_mg_dl, ldl_mg_dl, 
         tnf_a_ui_ml, adiponectin_ug_ml, crp_mg_dl)

d_meds <- med_data 

# PHENOTYPE DEFINITION #########################################################
all_data <- d_bq %>%
  left_join(d_anthropo, by = c("id", "visit")) %>%
  left_join(d_vs, by = c("id", "visit")) %>%
  left_join(d_meds, by = "id") %>%
    filter(antibiotics == 0) %>%
    select(-antibiotics) %>% 
  left_join(d_vol, by = "id") %>%
  left_join(stool_samples, by = c("id", "visit")) %>%
    filter(sample == 1) %>%
    select(-sample) %>%
  # Define obesity / non-obesity labels
  mutate(obese = if_else(bmi_kg_m2 >= 30, "O", "NO")) %>%
  # Define alterations in markers
  mutate(blood_pressure = 
           if_else(
             systolic_blood_pressure_mmhg > 130 |
               dyastolic_blood_pressure_mmhg > 85 |
               antihypertensives.x == 1 |
               antihypertensives.y == 1 |
               hypertension == 1,
             1, 0),
         hdl_level = if_else(
           sex == "Female" & hdl_mg_dl <= 50 |
             sex == "Male" & hdl_mg_dl <= 40,
           1, 0),
         trig_level = if_else(
           tri_mg_dl > 150 |
             triglyceridemia == 1 |
             cholesterol_lowering == 1 |
             statins == 1 |
             hypercholesterolemia == 1,
           1, 0),
         blood_glucose = if_else(
           glu_mg_dl > 100 |
             antidiabetes == 1 |
             diabetes == 1,
           1, 0)) %>%
  mutate(number_markers = 
           blood_pressure + blood_glucose + hdl_level + trig_level) %>%
  mutate(met_health = if_else(number_markers > 0, "MU", "MH")) %>%
  mutate("MetObesity" = paste0(met_health, obese)) %>%
  filter(!is.na(met_health))
summary(as.factor(all_data$MetObesity))

vols.mh.dup <- all_data %>%
  select(id, MetObesity) %>%
  distinct() %>%
  group_by(id) %>%
  summarise(n = n()) %>%
  filter(n > 1) %>% 
  pull(id)

df.dup <- all_data %>% 
  filter(id %in% vols.mh.dup) %>%
  select(id, visit, MetObesity) %>% 
  spread(id, MetObesity) %>% 
  column_to_rownames(var = "visit") %>% 
  t() %>% 
  data.frame() %>% 
  mutate(visit = case_when(
    X1 == "MHO" ~  1,
    X3 == "MHO" ~ 3,
    X1 == "MHNO" ~ 1,
    X3 == "MHNO" ~ 3,
    X1 == "MUO" ~ 3)) %>%
  rownames_to_column(var = "id") %>%
  select(id, visit)

df.part.1 <- all_data %>%
  right_join(df.dup, by = c("id", "visit"))

vols.one.visit <- all_data %>% 
  filter(!id %in% vols.mh.dup) %>% 
  select(id, visit) %>% 
  group_by(id) %>% 
  summarise(n = n()) %>% 
  filter(n == 1) %>% 
  pull(id)

df.part.2 <- all_data %>% filter(id %in% vols.one.visit)

df.part.3 <-  all_data %>% 
  filter(!id %in% vols.mh.dup) %>%
  filter(!id %in% vols.one.visit) %>%
  filter(visit == 1)


df.final <- bind_rows(df.part.1, df.part.2, df.part.3)
# sanity check: 1 muestra x participant:
length(unique(df.final$id)) == length(df.final$id) # OK!
summary(as.factor(df.final$MetObesity))

write_csv(df.final, "./ai4food_metadata_20250821.csv")