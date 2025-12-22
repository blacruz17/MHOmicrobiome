library(curatedMetagenomicData)
library(tidyverse)
library(janitor)
library(readxl)
library(gtsummary)

rm(list = ls())

# selected IDs from curatedMetagenomicsData:
id.df <- read_csv("SuppTable_cmd.csv") %>%
  filter(study_name == "MetaCardis_2020_a") %>%
  select(-study_name)

# Supplementary Tables from source publication:
st1a <- read_excel("41586_2021_4177_MOESM3_ESM.xlsx",
                   sheet = "ST 1a") %>%
  clean_names()

st1b <- read_excel("41586_2021_4177_MOESM3_ESM.xlsx",
                   sheet = "ST 1b") %>%
  clean_names() %>%
  rename(SampleID = sample_id)

# Supplementary material from Zenodo repository:
ze.demo  <- read.table("./metadata_zenodo/demographic_20201210.r",
                       header = TRUE)
ze.antb  <- read.table("./metadata_zenodo/antibiotics_20201210.r",
                       header = TRUE)
ze.drugs <- read.table("./metadata_zenodo/cmd_drugs_20201210.r",
                       header = TRUE)


# DATA WRANGLING ###############################################################
paper.data <- st1b %>%
  left_join(ze.demo,  by = "SampleID") %>%
  left_join(ze.antb,  by = "SampleID") %>%
  left_join(ze.drugs, by = "SampleID") %>%
  mutate(subject_id = paste0("M0", SampleID)) %>%
  inner_join(id.df, by = "subject_id")


all(paper.data$ANTIBIOTICS_TOTAL == 0)
# perfect!!
paper.data.2 <- paper.data %>% select(-colnames(ze.antb)[2:23])

count(paper.data, center)

summary(paper.data$age)

paper.data.MHNO <- paper.data.2 %>%
  filter(MetObesity == "MHNO") 

summary(as.factor(paper.data.MHNO$patient_group))

st1a[1, 10]
# Control group, subjects exempt of CAD, T2D, MS and with BMI<25


summary(as.factor(paper.data.MHNO$DRUGTOTAL))
paper.data.MHNO.nc <- paper.data.MHNO %>% remove_constant() 
colnames(paper.data.MHNO)[!colnames(paper.data.MHNO) %in% colnames(paper.data.MHNO.nc) ]
colnames(paper.data.MHNO.nc)

# METADATA INPUT_FEATURES ######################################################
input.features <- read.table("./input_features_zenodo/hub.lipo.v3.data.frame.r",
                             header = TRUE,
                             sep = "\t") %>%
  filter(SampleID %in% paper.data$SampleID)

input.features.2 <- input.features %>%
  filter(Feature %in% c("TPTG", # total plasma triglycerides mg/dL
                        "HDCH", # hdl cholesterol mg/dL
                        "LDCH", # ldl cholesterol mg/dL
                        "TPCH"  # total plasma cholesterol mg/dL
                        )) %>%
  select(-FeatureDisplayName) %>%
  pivot_wider(names_from = Feature,
              values_from = FeatureValue)

paper.data <- paper.data %>%
  left_join(input.features.2, by = "SampleID")

paper.data %>%
  select(patient_group, gender, bmi, age, center, 
         HDCH, TPTG,
         ANTIBIOTICS_TOTAL,
         GENDER, BMI_C, AGE, CENTER_C,
         ANTILIPID_C, ANTIHTA_C, ANTIDIAB_C, DRUGTOTAL,
         MetObesity) %>%
  tbl_summary(by = "MetObesity",
              type = list(gender = "dichotomous"))


paper.data %>% filter(MetObesity == "MHNO") %>% pull(TPTG) %>% sort()

paper.data %>%
  filter(MetObesity == "MHNO") %>% 
  select(HDCH, gender) %>% 
  mutate(hdl_ok = as.factor(
           case_when(gender == 0 & HDCH > 50 ~ "OK",
                     gender == 1 & HDCH > 40 ~ "OK",
                     is.na(HDCH) ~ NA,
                     .default = "MU"))) %>%
  pull(hdl_ok) %>%
  summary()

paper.data <- paper.data %>%
  mutate(hdl_ok = as.factor(
    case_when(gender == 0 & HDCH > 50 ~ "OK",
              gender == 1 & HDCH > 40 ~ "OK",
              is.na(HDCH) ~ NA,
              .default = "MU")),
    tg_ok = as.factor(
      case_when(TPTG <= 150 ~ "OK",
                is.na(TPTG) ~ NA,
                .default = "MU"))) %>%
  mutate(MU = as.factor(
    if_else(hdl_ok == "OK" & tg_ok == "OK", "OK", "MU"))) %>%
  mutate(MetObesity  = if_else(MetObesity == "MHNO" & MU == "MU", 
                               "MUNO",
                               MetObesity)) %>%
  select(-hdl_ok, -tg_ok, -MU) %>%
  rename(hdl = HDCH,
         ldl = LDCH,
         triglycerides = TPTG,
         cholesterol = TPCH)


# EXPORT TABLE #################################################################
write_csv(paper.data, "./metacardis_metadata_20250820.csv")

