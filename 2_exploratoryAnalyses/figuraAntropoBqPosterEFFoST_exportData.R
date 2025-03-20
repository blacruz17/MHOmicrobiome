rm(list = ls())

library(curatedMetagenomicData)
library(tidyverse)
library(janitor)
library(readxl)
library(ComplexHeatmap)
library(gtsummary)
library(RColorBrewer)
library(phyloseq)
library(rstatix)
library(ggpubr)
library(gt)



karlsson <- read_excel("suppKarlsson_PMID23719380.xlsx",
                   sheet = "Supplementary Table 3",
                   skip = 1) %>%
  clean_names()
karlsson <- karlsson[1:145, 1:30]  

karlsson <- karlsson %>%
  mutate(obese = if_else(bmi_kg_m2 >= 30, 'O', 'NO', 'missing'),
         met_health = if_else(
           classification == 'NGT' & # normal glucose transport
             insulin_y_n == "N" &
             statins_y_n == "N" &
             (oral_anti_diabetic_medication_no_medication_met_metformin_sulph_sulphonylurea == "-" |
                is.na(oral_anti_diabetic_medication_no_medication_met_metformin_sulph_sulphonylurea)) &
             fasting_glucose_mmol_l <= 6.1 &
             triglycerides_mmol_l <= 1.7 &
             hdl_mmol_l > 1.3,
           'MH', 'MU', 'missing'),
         subject_id = paste0("S", sample_id),
         sample_id = subject_id)  %>%
  mutate("MetObesity" = as.factor(paste0(met_health, obese)))

summary(karlsson$MetObesity)

feng <- sampleMetadata %>%
  filter(study_name == "FengQ_2015")  %>%
  filter(study_condition == "control")  %>% 
  clean_names() %>%
  select(where(~ !all(is.na(.x))))

feng <- read_excel("suppFeng_PMID25758642.xlsx",
                   sheet = "SD 1",
                   skip = 3,
                   na = c("", "NA", "n.a.")) %>%
          clean_names()

feng <- feng %>%
  filter(state == "controls") %>% 
  mutate(obese = if_else(bmi >= 30, 'O', 'NO', 'missing')) %>%
  mutate(met_health = if_else(
    # medicacion / enfermedad
    fatty_liver_in_ultrasound_1_yes_0_no == 0 &
      diabetes_1_yes_0_no == 0 &
      hypertension_1_yes_0_no == 0 &
      met_s_1_yes_0_no == 0 &
      # parametros bioquimicos
      fasting_glucose_mg_l <= 100 &
      tg_mg_l <= 150 &
      ((gender_1_male_2_female == 1 & hdl_mg_l > 40)|
         (gender_1_male_2_female == 2 & hdl_mg_l > 50)),
    'MH', 'MU', 'missing'))  %>%
  mutate("MetObesity" = as.factor(paste0(met_health, obese))) %>%
  mutate(sample_id = paste0("SID", x1)) %>%
  mutate(subject_id = sample_id)

summary(feng$MetObesity)

hmp <- sampleMetadata %>%
  filter(study_name == "HMP_2019_t2d",
         is.na(antibiotics_current_use))  %>% 
  clean_names() %>%
  select(where(~ !all(is.na(.x)))) %>%
  mutate(met_health = if_else(
    # medicacion / enfermedad
    disease == "healthy" &
      # parametros bioquimicos
      glucose <= 100 &
      triglycerides <= 150 &
      ((gender == "male" & hdl > 40)|
         (gender == "female" & hdl > 50)),
    'MH', 'MU', 'missing')) %>%
  filter(met_health != "missing",
         !disease %in% c("IGT;respiratoryinf", 
                         "respiratoryinf",
                         "T2D;respiratoryinf"))

hmp  %>% 
  filter(met_health != "missing") %>% 
  select(subject_id, met_health) %>% 
  distinct() %>% 
  group_by(subject_id) %>% 
  summarise(n = n()) %>% 
  filter(n>1)

hmp.rep <- hmp  %>% 
  filter(met_health != "missing") %>% 
  select(subject_id, met_health) %>% 
  distinct() %>% 
  group_by(subject_id) %>% 
  summarise(n = n()) %>% 
  filter(n>1) %>%
  pull(subject_id)

hmp.part1 <- hmp %>% filter(subject_id %in% hmp.rep)
hmp.part1 <- hmp.part1 %>% 
  filter(sample_id == "HMP2_J00883_M_ST_T0_B0_0120_ZM7JY3G-02_HA986ADXX")

hmp.part2 <- hmp %>% 
  filter(!subject_id %in% hmp.rep) %>% 
  distinct(subject_id, .keep_all = TRUE)

hmp <- rbind(hmp.part1, hmp.part2)

summary(as.factor(hmp$subject_id))

hmp_extra <- read_excel("suppHMP_PMID31142858.xlsx",
                        sheet = "S1_Subjects") %>%
  clean_names() %>%
  mutate(obese = if_else(bmi >= 30, "O", "NO"))

hmp <- hmp %>% 
  left_join(hmp_extra, by = "subject_id") %>%
  mutate(MetObesity = as.factor(paste0(met_health, obese)))
summary(hmp$MetObesity)

metacard <- sampleMetadata %>%
  filter(study_name == "MetaCardis_2020_a") %>%
  clean_names() %>%
  select(where(~ !all(is.na(.x)))) %>%
  mutate_if(is.character, as.factor) %>% 
  filter(antibiotics_current_use == "no") %>% 
  mutate(obese = if_else(bmi > 30, "O", "NO"),
         met_health = if_else(disease == "healthy" & 
                                treatment == "no" & 
                                triglycerides <= 150, 
                              "MH", "MU"),
         MetObesity = as.factor(paste0(met_health, obese))) %>%
  filter(MetObesity %in% c("MHO", "MHNO", "MUO", "MUNO"))
summary(as.factor(metacard$MetObesity))

metacard %>%
    separate_rows(treatment, sep = ";") %>% pull(treatment) %>% unique()
metacard %>% select(study_condition, disease) %>% distinct()

ages.mc <- readxl::read_excel("suppMetaCardis_PMID35177860_limpio.xlsx") %>%
  
  rename(sample_id = id, age = age_years) %>%
  select(sample_id, age) %>% # "x12MCx3370" 
  mutate(sample_id = paste0("M0", sample_id)) %>%
  distinct()

metacard <- metacard %>% left_join(ages.mc, by = "sample_id")


## Karlsson --------------------------------------------------------------------
karlsson <- karlsson %>%
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
feng <- feng %>% 
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
  # crear columna disease
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
hmp <- hmp[ , c(1:3, 5:7, 9, 21, 23:48, 50:57)] %>%
  select(-c(study_condition, location, population, consented, class, ethnicity, 
            adj_age, num_all_visits, days_span, num_healthy_visits, sspg)) %>%
  rename(gender = gender.x,
         albumin = albumine)

## MetaCardis ------------------------------------------------------------------
metacardis <- metacard %>%
  select(-c(study_name, 
            body_site, 
            antibiotics_current_use,
            study_condition, age_category, non_westernized,
            sequencing_platform, pmid, ncbi_accession, curator, location,
            disease_subtype, smoke, number_bases)) %>%
  rename(ldl = ldl_2,
         hscrp = hs_crp) %>%
  filter(MetObesity %in% c("MHO", "MHNO", "MUO", "MUNO"))

ai4food <- read_csv("tabla_mho_ai4food.csv") %>%
  rename(subject_id = id_voluntario,
         gender = sexo) %>%
  mutate(sample_id = subject_id,
         subject_id = paste0(as.character(subject_id), "_A4F"),
         country = "ESP") %>%
  select(-c(glu_mg_dl, hdl_mg_dl, tri_mg_dl, 
            blood_pressure, hdl_level, trig_level, blood_glucose,
            number_markers)) 

ai4food_2 <- read_tsv("../datos/tablas_export_v2sql/bioquimicas_v2.tsv",
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

ai4food_3 <- read_tsv("../datos/tablas_export_v2sql/datos_antropometricos_v2.tsv",
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

ai4food$treatment <- "no"
ai4food <- ai4food %>%
  mutate(treatment = if_else(hipotensores == 1, 
                             "antihypertensive", 
                             treatment)) %>%
  mutate(treatment = if_else(antidiabetico == 1, 
                             "antidiabetic",
                             treatment)) %>%
  mutate(treatment = if_else(estatina == 1, 
                             if_else(treatment != "no",
                                     paste0(treatment, ";statins"),
                                     "statins"),
                             treatment)) %>%
  mutate(treatment = if_else(reductor_de_colesterol == 1, 
                             if_else(treatment != "no",
                                     paste0(treatment, ";cholesterol"),
                                     "cholesterol"),
                             treatment))


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

ages.ai4food <- readxl::read_excel("../datos/brutos_edad_limpio.xlsx") %>%
  
  rename(subject_id = vol_id, age = muestra) %>%
  select(subject_id, age) %>%
  mutate(subject_id = paste0(as.character(subject_id), "_A4F"))

ai4food <- ai4food %>% left_join(ages.ai4food, by = "subject_id")

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
                           region = "Europe",
                           bmi = as.numeric(bmi)), 
                  metacardis %>%
                    mutate(study_name = "MetaCardis_2020_a",
                           region = "Europe"), 
                  ai4food %>%
                    mutate(study_name = "AI4Food",
                           region = "Europe",
                           sample_id = as.character(sample_id)))
char_cols <- sapply(data, is.character)
data[char_cols][data[char_cols] == '-'] <- NA

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


data <- data %>%
  mutate(gender = case_when(gender == "H" ~ "male",
                            gender == "M" ~ "female",
                            .default = gender),
         disease = if_else(disease %in% c("control", "no"), "healthy", disease))

og_treatment <- data$treatment == "no"
data$treatment <- str_count(data$treatment, ";") + 1
data$treatment[og_treatment] <- 0


data <- data %>%
  mutate(
    hip = if_else(hip < 50, NA, hip),
    whr = if_else(whr > 5, NA, whr)
  ) %>%
  mutate(tnfa_1 = if_else(study_name == "AI4Food", tnfa, NA),
         tnfa_2 = if_else(study_name == "KarlssonFH_2013", tnfa, NA)) %>%
  select(-tnfa)


## ----subsetSamples---------------------------------------------------------------------------------------------------
physeq <- readRDS("adj_physeq_MMUPHIN_mpa30_20240911.rds")
data <- data %>% filter(sample_id %in% sample_names(physeq))
rm(physeq)


## ----gtsummary-------------------------------------------------------------------------------------------------------

tbl_todo <- data %>%
  filter(!is.na(met_health), !is.na(obese),
         MetObesity != "missingmissing",
         !is.na(gender) # estorba un poco
  ) %>%
  mutate(treatment = as.factor(if_else(treatment == 0, "No", "Yes"))) %>%
  select(MetObesity,
         study_name,
         gender, age, treatment, bmi, 
         weight, hip, waist, whr,
         glucose, hba1c, homa, insulin,
         triglycerides, hdl, ldl, cholesterol, 
         diastolic, systolic,
         tnfa_1, tnfa_2,
         adiponectin, albumin, crp) %>%
  tbl_summary(by = "MetObesity",
              missing = "no") %>%
  add_n() %>%
  add_p() %>%
  add_q() %>%
  bold_p(q = TRUE)
tbl_todo


## ----pairwise--------------------------------------------------------------------------------------------------------
dunn.tb <- data %>%
  filter(!is.na(met_health), !is.na(obese),
         MetObesity != "missingmissing",
         !is.na(gender) # estorba un poco
  ) %>%
  select(MetObesity,
         age, bmi, 
         weight, hip, waist, whr,
         glucose, hba1c, homa, insulin,
         triglycerides, hdl, 
         diastolic, systolic,
         tnfa_1,
         adiponectin, crp) %>% 
  reshape2::melt(id.vars = "MetObesity") %>% 
  group_by(variable) %>%
  dunn_test(value ~ MetObesity,
            p.adj = "BH")

dunn.tb

## ----wranglingPvals--------------------------------------------------------------------------------------------------
stat.test <- data %>%
  filter(!is.na(met_health), !is.na(obese),
         MetObesity != "missingmissing",
         !is.na(gender) # estorba un poco
  ) %>%
  select(MetObesity,
         age, bmi, 
         weight, hip, waist, whr,
         glucose, hba1c, homa, insulin,
         triglycerides, hdl, 
         diastolic, systolic,
         tnfa_1,
         adiponectin, crp) %>% 
  reshape2::melt(id.vars = "MetObesity") %>% 
  mutate(variable2 = variable) %>%
  group_by(variable2) %>%
  dunn_test(value ~ MetObesity,
            p.adj = "BH") %>%
  add_xy_position(scales = "free_y")


pal <- c("#264653", "#2A9D8F", "#703d57", "#edafb8", "#e6af2e")

var.labs <- c("Age", "BMI", 
              "Weight", "Hip", "Waist", "Waist-hip ratio",
              "Glucose", "Hb1Ac", "HOMA-IR", "Insulin",
              "Triglycerides", "HDL", 
              "Diastolic blood pressure", "Systolic blood pressure",
              "TNF- \U03B1 (AI4Food)", 
              "Adiponectin", "C-reactive protein")
names(var.labs) <- c("age", "bmi", 
                     "weight", "hip", "waist", "whr",
                     "glucose", "hba1c", "homa", "insulin",
                     "triglycerides", "hdl", 
                     "diastolic", "systolic",
                     "tnfa_1",
                     "adiponectin", "crp")



## ----statPlot--------------------------------------------------------------------------------------------------------

p.metadata_exploration <- data %>%
  filter(!is.na(met_health), !is.na(obese),
         MetObesity != "missingmissing",
         !is.na(gender)) %>% # estorba un poco
  mutate(treatment = as.factor(if_else(treatment == 0, "No", "Yes"))) %>%
  select(MetObesity,
         age, bmi, 
         weight, hip, waist, whr,
         glucose, hba1c, homa, insulin,
         triglycerides, hdl, 
         diastolic, systolic,
         tnfa_1,
         adiponectin, crp) %>%
  reshape2::melt(id.vars = "MetObesity") %>%
  rename(variable2 = variable) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = MetObesity, y = value)) +
  geom_boxplot(aes(fill = MetObesity), 
               alpha=.6,
               outlier.shape = NA,
               outlier.size = 2,
               color = "black",
               # width = .1
  ) +
  geom_jitter(aes(shape = MetObesity, color = MetObesity),
              width = .2,
              size = .8) +
  scale_color_manual(values = pal,
                     aesthetics = c("fill", "color"),
                     name = "") +
  facet_wrap(. ~ variable2, 
             ncol = 3,
             scales = "free_y",
             labeller = labeller(variable2 = var.labs)) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  stat_pvalue_manual(stat.test,
                     hide.ns = TRUE,
                     tip.length = 0.01,
                     bracket.shorten = 0.05,
                     # step.increase = .1,
                     # bracket.size = 0.3,
                     label.size = 2.5
  ) +
  theme_linedraw() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    plot.title = element_text(face = 'bold'),
    axis.text = element_text(size = 6),
    strip.text = element_text(face = 'bold', size = 7,
                              color = "black"),
    strip.background = element_rect(
      color="white", fill="white",  linetype="solid"
    ))
p.metadata_exploration


## ----saveFig---------------------------------------------------------------------------------------------------------
ggsave("./metadataExploration_DunnTest_20250225.png",
       plot = p.metadata_exploration,
       bg = "white",
       height = 20,
       width = 14,
       units = "cm",
       dpi = 1200)


## ----otherAnthropoPars-----------------------------------------------------------------------------------------------
data %>%
  filter(!is.na(met_health), !is.na(obese),
         MetObesity != "missingmissing",
         !is.na(gender)) %>% 
  select(MetObesity, bmi) %>%
  ggplot(aes(x = MetObesity, y = bmi)) +
  geom_boxplot(aes(fill = MetObesity), 
               alpha=.6,
               outlier.shape = NA,
               outlier.size = 2,
               color = "black"
  ) +
  geom_jitter(aes(shape = MetObesity, color = MetObesity),
              width = .2,
              size = .8) +
  scale_color_manual(values = pal,
                     aesthetics = c("fill", "color"),
                     name = "") +
  theme_bw() +
  theme(
    legend.position = "none"
  )


data %>%
  filter(!is.na(met_health), !is.na(obese),
         MetObesity != "missingmissing",
         !is.na(gender)) %>% 
  select(MetObesity, waist, hip, whr) %>%
  reshape2::melt(id.vars = "MetObesity") %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = MetObesity, y = value)) +
  geom_boxplot(aes(fill = MetObesity), 
               alpha=.6,
               outlier.shape = NA,
               outlier.size = 2,
               color = "black"
  ) +
  geom_jitter(aes(shape = MetObesity, color = MetObesity),
              width = .2,
              size = .8) +
  scale_color_manual(values = pal,
                     aesthetics = c("fill", "color"),
                     name = "") +
  theme_bw() +
  theme(
    legend.position = "none"
  ) +
  facet_wrap(. ~variable, scales = "free_y")

data %>%
  filter(!is.na(met_health), !is.na(obese),
         MetObesity != "missingmissing",
         !is.na(gender)) %>%
  select(MetObesity, bmi) %>%
  mutate(class = case_when(bmi <= 18.5            ~ "Underweight",
                           bmi >= 18.5 & bmi < 25 ~ "Normal weight",
                           bmi >= 25   & bmi < 30 ~ "Overweight",
                           bmi >= 30   & bmi < 35 ~ "ClassI",
                           bmi >= 35   & bmi < 40 ~ "ClassII",
                           bmi >= 40              ~ "ClassIII")) %>%
  mutate(class = as.factor(class)) %>%
  group_by(MetObesity) %>%
  mutate(total = n()) %>%
  group_by(MetObesity, class) %>%
  mutate(classGroup = n()) %>%
  arrange(MetObesity) %>%
  ungroup() %>%
  mutate(perClassGroup = 100 * classGroup/total) %>%
  select(MetObesity, class, classGroup, perClassGroup) %>%
  distinct()

