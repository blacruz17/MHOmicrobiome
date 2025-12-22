library(tidyverse)
library(janitor)
library(readxl)

rm(list = ls())


setwd(iniDir)

# KARLSSON FH 2013 ---------
# Supplementary Table from source publication:
karlsson <- read_excel(paste0(dataDir, "suppKarlsson_PMID23719380.xlsx"),
                       sheet = "Supplementary Table 3",
                       skip = 1) %>%
  clean_names()
karlsson <- karlsson[1:145, 1:30] 

karlsson <- karlsson %>%
  mutate(obese = if_else(bmi_kg_m2 >= 30, 'O', 'NO', 'missing'),
         met_health = if_else(
           classification == 'NGT' &
             insulin_y_n == "N" &
             statins_y_n == "N" &
             (oral_anti_diabetic_medication_no_medication_met_metformin_sulph_sulphonylurea == "-" |
                is.na(oral_anti_diabetic_medication_no_medication_met_metformin_sulph_sulphonylurea)) &

             fasting_glucose_mmol_l <= 6.1 &
             triglycerides_mmol_l <= 1.7 &
             hdl_mmol_l > 1.3, # only women in this study
           'MH', 'MU', 'missing'),
         subject_id = paste0("S", sample_id),
         sample_id = subject_id)  %>%
  mutate(hip_cm = wc_cm/whr_cm_cm) %>%
  mutate("MetObesity" = as.factor(paste0(met_health, obese)))

summary(karlsson$MetObesity)

write_csv(karlsson, "karlsson_metadata_20250820.csv")

# FENG Q 2015 ----
# Supplementary Table from source publication:
feng <- read_excel(paste0(dataDir, "suppFeng_PMID25758642.xlsx"),
                   sheet = "SD 1",
                   skip = 3,
                   na = c("", "NA", "n.a.")) %>%
  clean_names()

feng <- feng %>%
  filter(state == "controls") %>% # IMPORTANT, otherwise adenoma or crc!
  mutate(obese = if_else(bmi >= 30, 'O', 'NO', 'missing')) %>%
  mutate(met_health = if_else(

    fatty_liver_in_ultrasound_1_yes_0_no == 0 &
      diabetes_1_yes_0_no == 0 &
      hypertension_1_yes_0_no == 0 &
      met_s_1_yes_0_no == 0 &

      fasting_glucose_mg_l <= 100 &
      tg_mg_l <= 150 &
      ((gender_1_male_2_female == 1 & hdl_mg_l > 40)|
         (gender_1_male_2_female == 2 & hdl_mg_l > 50)),
    'MH', 'MU', 'missing'))  %>%
  mutate("MetObesity" = as.factor(paste0(met_health, obese))) %>%
  mutate(sample_id = paste0("SID", x1)) %>%
  mutate(subject_id = sample_id)

summary(feng$MetObesity)

write_csv(feng, "feng_metadata_20250820.csv")
