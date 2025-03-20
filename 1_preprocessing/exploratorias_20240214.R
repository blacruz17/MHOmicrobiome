library(tidyverse)
library(janitor)
library(reshape2)
library(readxl)
library(RColorBrewer)
library(ComplexHeatmap)

# DATA IMPORT  ###############################################################
## Karlsson --------------------------------------------------------------------
karlsson <- read_tsv("karlsson_metadata.tsv") %>%
  rename(oral = oral_anti_diabetic_medication_no_medication_met_metformin_sulph_sulphonylurea) %>%
  mutate(treatment = case_when(
    insulin_y_n == "Y" & statins_y_n == "Y" ~ "insulin;statins",
    insulin_y_n == "Y" & (statins_y_n == "N" | is.na(statins_y_n)) ~ "insulin",
    statins_y_n == "Y" & (insulin_y_n == "N" | is.na(insulin_y_n)) ~ "statins")) %>%
  mutate(treatment = case_when(is.na(treatment) & is.na(oral) ~ NA,
                               is.na(treatment) & oral == "-" ~ "no",
                               is.na(treatment) & oral != "-" ~ oral,
                               !is.na(treatment) & oral == "-" ~ treatment,
                               !is.na(treatment) & oral != "-" ~
                                 paste(treatment, oral, sep = ";")),
         gender = "female") %>%
  rename(waist = wc_cm,
         ygt = y_gt_mkat_l,
         glucose = fasting_glucose_mmol_l,
         insulin = fasting_insulin_m_u_l,
         hba1c = hb_a1c_mmol_mol,
         hscrp = hs_crp_mg_l,
         cpeptide = c_peptide_nmol_l,
         tnfa = tn_fa_ng_l,
         il1 = il_1_pg_ml,
         age = age_years,
         gad = gad_antibodies_units1,
         bmi = bmi_kg_m2,
         whr = whr_cm_cm,
         cholesterol = cholesterol_mmol_l,
         triglycerides = triglycerides_mmol_l,
         hdl = hdl_mmol_l,
         ldl = ldl_mmol_l,
         creatinine = creatinine_mmol_l,
         adiponectin = adiponectin_mg_l,
         lepin = leptin_mg_l,
         glp = glp_1_pmol_l,
         fgf = fgf_19_pg_ml) %>%
  select(-c(country_of_birth, oral, insulin_y_n, statins_y_n,
            years_in_sweden)) %>%
  rename(disease = classification) %>%
  mutate(glucose = as.numeric(glucose))

## Feng ------------------------------------------------------------------------
feng <- read_tsv("feng_metadata.tsv",
                 na = c("", "NA", "n.a")) %>% 
  rename(age = age_yrs,
         waist = waist_cm,
         hip = hip_cm,
         whr = whr_waist_hip_ratio,
         gender = gender_1_male_2_female,
         fattyliver = fatty_liver_in_ultrasound_1_yes_0_no,
         ggt = ggt_u_l,
         got = got_ast_u_l,
         gpt = gpt_alt_u_l,
         insulin = fasting_insulin_m_u_m_l,
         glucose = fasting_glucose_mg_l,
         homa = homa_index,
         diabetes = diabetes_1_yes_0_no,
         hba1c = hba1c_percent,
         crp = crp_mg_l,
         ferritin = ferritin_ng_m_l,
         hb = hb_g_l,
         triglycerides = tg_mg_l,
         hdl = hdl_mg_l,
         ldl = ldl_mg_l,
         hypertension = hypertension_1_yes_0_no,
         mets = met_s_1_yes_0_no) %>%
  mutate(subject_id = as.character(subject_id),
         whr = as.numeric(whr),
         gender = if_else(gender == 1, "male", "female")) %>%
  # create disease column
  mutate(disease = case_when(
    diabetes == 0 & hypertension == 0 & mets == 0 & fattyliver == 1 ~ "fattyliver",
    diabetes == 0 & hypertension == 1 & mets == 1 & fattyliver == 0 ~ "hypertension;MetS",
    diabetes == 0 & hypertension == 1 & mets == 1 & fattyliver == 1 ~ "hypertension;MetS;fattyliver",
    diabetes == 0 & hypertension == 0 & mets == 0 & fattyliver == 0 ~ "no",
    diabetes == 1 & hypertension == 1 & mets == 1 & fattyliver == 1 ~ "diabetes;hypertension;MetS;fattyliver",
    diabetes == 0 & hypertension == 1 & mets == 0 & fattyliver == 1 ~ "hypertension;fattyliver",
    diabetes == 0 & hypertension == 1 & mets == 0 & fattyliver == 0 ~ "hypertension",
    diabetes == 1 & hypertension == 0 & mets == 0 & fattyliver == 1 ~ "diabetes;fattyliver",
    diabetes == 1 & hypertension == 1 & mets == 1 & fattyliver == 0 ~ "diabetes;hypertension;MetS",
    diabetes == 1 & hypertension == 0 & mets == 1 & fattyliver == 1 ~ "diabetes;MetS;fattyliver")) %>%
  select(-c(diabetes, hypertension, mets, fattyliver))

## HMP 2019 --------------------------------------------------------------------
hmp <- read_tsv("hmp_metadata.tsv")
hmp <- hmp[ , c(1:3, 5:7, 9, 21, 23:48, 50:57)] %>%
  select(-c(study_condition, location, population, consented, class, ethnicity, 
            adj_age, num_all_visits, days_span, num_healthy_visits, sspg)) %>%
  rename(gender = gender.x,
         albumin = albumine)

## MetaCardis ------------------------------------------------------------------
metacardis <- read_tsv("metacardis_metadata.tsv") %>%
  select(-c(study_name, body_site, antibiotics_current_use,
            study_condition, age_category, non_westernized,
            sequencing_platform, pmid, ncbi_accession, curator, location,
            disease_subtype, smoke, number_bases)) %>%
  rename(ldl = ldl_2,
         hscrp = hs_crp) %>%
  filter(MetObesity %in% c("MHO", "MHNO", "MUO", "MUNO"))

## JieZ ------------------------------------------------------------------------
jiez <- read_tsv("jiez_metadata.tsv") %>%
  rename(disease = coronary_type,
         bmi = body_mass_index_bmi,
         whr = waist_hip_ratio_whr,
         triglycerides = trig_mmol_l,
         ldl = ldlc_mmol_l,
         cholesterol = chol_mmol_l,
         hdl = hdlc_mmol_l,
         hba1c = hb_a1c_percent,
         albumin = alb_u_l,
         creatinine = crea_umol_l,
         glucose = fbg) %>%
  select(-smoking)
jiez <- jiez[ ,c(1, 3:36, 56, 68:71)] 
colnames(jiez)[3:35] <- gsub("_.*", "", colnames(jiez)[3:35])

## QinJ ------------------------------------------------------------------------
qinj <- read_tsv("qinj_metadata.tsv") %>%
  select(-c(body_site, antibiotics_current_use,
            study_condition, age_category, sbp_mm_hg, dbp_mm_hg,
            diabetic_y_or_n, fcp_ng_ml, tcho_mmol_l, bmi_kg_m2)) %>%
  rename(cpeptide = c_peptide,
         glucose = fbg_mmol_l,
         hba1c = hbalc_percent,
         triglycerides = tg_mmol_l,
         insulin = fins_m_u_l,
         hdl = hdl_mmol_l,
         ldl = ldl_mmol_l,
         diastolic = dyastolic_p,
         height = height_cm,
         systolic = systolic_p,
         weight = weight_kg)

qinj <- qinj[, c(1:5, 16:33)]


## YuJ_2015 -------------------------------------------------------
yuj <- read_tsv("yuj_metadata.tsv") 
yuj <- yuj[, c(1:3, 7, 8, 10, 22, 24:26, 29:31, 34:36)] %>%
  rename(glucose = fasting_glucose)


## AI4Food ---------------------------------------------------------------------
ai4food <- read_csv("tabla_mho_ai4food.csv") %>%
  rename(subject_id = id_voluntario,
         gender = sexo) %>%
  mutate(sample_id = subject_id,
         subject_id = paste0(as.character(subject_id), "_A4F"),
         country = "ESP") %>%
  select(-c(glu_mg_dl, hdl_mg_dl, tri_mg_dl, 
            blood_pressure, hdl_level, trig_level, blood_glucose,
            number_markers)) 

ai4food_2 <- read_tsv("./biochem_v2.tsv",
                      na = c("", "NA", "\\N")) %>%
  select(id_voluntario, visita,
         hba1c_perc, glu_mg_dl, alb_g_dl, chol_mg_dl, tri_mg_dl,
         hdl_mg_dl, ldl_mg_dl, pcr_mg_dl, insulina_uui_ml, homa,
         tnf_a_ui_ml, adiponectina_ug_ml) %>%
  rename(subject_id = id_voluntario,
         hba1c = hba1c_perc,
         glucose = glu_mg_dl,
         albumin = alb_g_dl,
         cholesterol = chol_mg_dl,
         triglycerides = tri_mg_dl,
         hdl = hdl_mg_dl,
         ldl = ldl_mg_dl,
         crp = pcr_mg_dl,
         insulin = insulina_uui_ml,
         tnfa = tnf_a_ui_ml,
         adiponectin = adiponectina_ug_ml) %>%
  mutate(subject_id = paste0(as.character(subject_id), "_A4F"),
         crp = as.numeric(gsub("<0.003", 0.003, crp)),
         hba1c = 100 * hba1c)

ai4food_3 <- read_tsv("./anthropo_v2.tsv",
                      na = c("", "NA", "\\N")) %>%
  mutate(visita = as.numeric(gsub(0, 1, visita))) %>%
  select(id_voluntario, visita,
         imc_kg_m2, presion_sistolica_mmhg, presion_diastolica_mmhg,
         cadera_cm, cintura_cm, peso_actual_kg) %>%
  rename(subject_id = id_voluntario,
         bmi = imc_kg_m2,
         systolic = presion_sistolica_mmhg,
         diastolic = presion_diastolica_mmhg,
         waist = cintura_cm,
         hip = cadera_cm,
         weight = peso_actual_kg) %>%
  mutate(subject_id = paste0(as.character(subject_id), "_A4F"),
         whr = waist/hip)

ai4food <- left_join(ai4food, ai4food_2,
                     by = c("subject_id", "visita")) %>%
            left_join(ai4food_3,
            by = c("subject_id", "visita"))

# create treatment column:
ai4food$treatment <- "no"
ai4food <- ai4food %>%
  mutate(treatment = if_else(hipotensores == 1, 
                             "antihypertensive", 
                             treatment)) %>%
  mutate(treatment = if_else(antidiabetico == 1, 
                             "antidiabetic",
                             treatment)) %>%
  # statins
  mutate(treatment = if_else(estatina == 1, 
                             if_else(treatment != "no",
                                     paste0(treatment, ";statins"),
                                     "statins"),
                             treatment)) %>%
  # cholesterol lowering
  mutate(treatment = if_else(reductor_de_colesterol == 1, 
                             if_else(treatment != "no",
                                     paste0(treatment, ";cholesterol"),
                                     "cholesterol"),
                             treatment))

# create disease column:
ai4food$disease <- "control"
ai4food <- ai4food %>%
  mutate(disease = if_else(hipercolesterolemia == 1,
                           "hypercholesterolemia",
                           disease)) %>%
  mutate(disease = if_else(trigliceridemia == 1,
                           "triglyceridemia",
                           disease)) %>%
  mutate(disease = if_else(diabetes == 1, 
                           "diabetes", 
                           disease))  %>%
  mutate(disease = if_else(hipertension == 1, 
                           "hypertension", 
                           disease)) 
ai4food <- ai4food[ , c(1, 14:40)]


# JOIN DATAFRAMES ##############################################################
data <- bind_rows(karlsson %>%
                    mutate(study_name = "KarlssonFH_2013",
                           adiponectin = as.numeric(gsub("-", NA, adiponectin)),
                           region = "Europe"), 
                  feng %>%
                    mutate(study_name = "FengQ_2015",
                           homa = as.numeric(homa),
                           crp = as.numeric(crp),
                           region = "Europe"), 
                  hmp %>%
                    mutate(study_name = "HMP_2019_t2d",
                           region = "Europe"), 
                  metacardis %>%
                    mutate(study_name = "MetaCardis_2020_a",
                           region = "Europe"), 
                  jiez %>%
                    mutate(study_name = "JieZ_2017",
                           region = "China"), 
                  qinj %>%
                    mutate(study_name = "QinJ_2012",
                           region = "China"), 
                  yuj %>%
                    mutate(study_name = "YuJ_2015",
                           region = "China"), 
                  ai4food %>%
                    mutate(study_name = "AI4Food",
                           region = "Europe",
                           sample_id = as.character(sample_id)))
# Identify columns with character values
char_cols <- sapply(data, is.character)

# Replace '-' with NA in character columns
data[char_cols][data[char_cols] == '-'] <- NA

## Unit changes #########################################################
check.nums <- c("cholesterol", "triglycerides", "hdl", "ldl", "glucose",
                "insulin", "hba1c", "cpeptide", "tnfa", "crp", "albumin")
p.before.units <- data[, check.nums] %>%
  cbind(data[, "study_name"]) %>%
  reshape2::melt(id.vars = "study_name") %>%
  group_by(study_name, variable) %>%
  filter(!is.na(value)) %>%
  summarise(mean = mean(value)) %>%
  arrange(variable) %>%
  ggplot(aes(x = study_name, y = mean, fill = study_name)) +
  geom_col(color = "black") +
  facet_wrap(. ~  variable) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
p.before.units

data <- data %>%
  mutate(cholesterol = if_else(study_name %in% c("KarlssonFH_2013", "JieZ_2017"),
                               cholesterol * 38.67, cholesterol),
         triglycerides = if_else(study_name %in% c("KarlssonFH_2013", "JieZ_2017", "QinJ_2012"),
                                 triglycerides * 88.57, triglycerides),
         hdl = if_else(study_name %in% c("KarlssonFH_2013", "JieZ_2017", "QinJ_2012"),
                       hdl * 38.67, hdl),
         ldl = if_else(study_name %in% c("KarlssonFH_2013", "JieZ_2017", "QinJ_2012"),
                       ldl * 38.67, ldl),
         glucose = if_else(study_name %in% c("KarlssonFH_2013", "JieZ_2017", "QinJ_2012"),
                           glucose * 18.018, glucose),
         hba1c = if_else(study_name == "KarlssonFH_2013", 
                         (hba1c * 0.0915) + 2.15, 
                         hba1c),
         albumin = if_else(study_name %in% c("HMP_2019_t2d", "JieZ_2017"),
                           albumin/10, albumin),
         crp = if_else(study_name == "JieZ_2017", crp/10, crp)
  )
# check:
p.after.units <- data[, check.nums] %>%
  cbind(data[, "study_name"]) %>%
  reshape2::melt(id.vars = "study_name") %>%
  group_by(study_name, variable) %>%
  filter(!is.na(value)) %>%
  summarise(mean = mean(value)) %>%
  arrange(variable) %>%
  ggplot(aes(x = study_name, y = mean, fill = study_name)) +
  geom_col(color = "black") +
  facet_wrap(. ~  variable) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), legend.position = "top")
p.after.units

## Categorical variables #######################################################
data <- data %>%
  mutate(gender = case_when(gender == "H" ~ "male",
                            gender == "M" ~ "female",
                            .default = gender),
         disease = if_else(disease %in% c("control", "no"), "healthy", disease))
og_treatment <- data$treatment == "no"
data$treatment <- str_count(data$treatment, ";") + 1
data$treatment[og_treatment] <- 0

## Location #####
data <- data %>%
  mutate(western = as.factor(if_else(study_name %in% c("QinJ_2012", 
                                                       "JieZ_2017",
                                                       "YuJ_2015"),
                                    "No", "Yes")),
         study_name = as.factor(study_name))

# Save ########################  
write_csv(data, "metadata_processed.csv")