## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
rm(list = ls())

library(tidyverse)
library(janitor)
library(readxl)
library(gtsummary)
library(RColorBrewer)
library(rstatix)
library(ggpubr)
library(gt)
library(emmeans)
library(phyloseq)


## ----ai4food-----------------------------------------------------------------------------------------------------------------------------------------
ai4food <- read_csv(
  "ai4food_metadata_20250821.csv"
  ) %>%
  select(id, sex, age, bmi_kg_m2,
         waist_cm, hip_cm,
         hba1c_perc, insulin_uui_ml, homa, glu_mg_dl, 
         chol_mg_dl, tri_mg_dl, hdl_mg_dl, ldl_mg_dl, 
         tnf_a_ui_ml, adiponectin_ug_ml, crp_mg_dl, bmi_kg_m2, 
         systolic_blood_pressure_mmhg, dyastolic_blood_pressure_mmhg,
         MetObesity) %>%
  mutate(whr = waist_cm/hip_cm)


## ----metaCardis--------------------------------------------------------------------------------------------------------------------------------------
metacardis <- read_csv(
  "metacardis_metadata_20250820.csv"
  )  %>%
  select(SampleID, gender, bmi, age, MetObesity, 
         ldl, hdl, triglycerides, cholesterol) %>%
  rename(id = SampleID,
         sex = gender,
         bmi_kg_m2 = bmi,
         chol_mg_dl = cholesterol, 
         tri_mg_dl = triglycerides, 
         hdl_mg_dl = hdl, 
         ldl_mg_dl = ldl) %>%
  mutate(sex = if_else(sex == 0, "Female", "Male"))

all(colnames(metacardis) %in% colnames(ai4food))


## ----Karlsson----------------------------------------------------------------------------------------------------------------------------------------
karlsson <- read_csv("karlsson_metadata_20250820.csv",
                     na = c("", "NA", "-")) %>%
  select(sample_id, age_years, bmi_kg_m2, 
         whr_cm_cm, wc_cm,
         cholesterol_mmol_l, triglycerides_mmol_l, hdl_mmol_l, ldl_mmol_l,
         fasting_glucose_mmol_l, fasting_insulin_m_u_l, hb_a1c_mmol_mol,
         adiponectin_mg_l, tn_fa_ng_l, MetObesity) %>%
  rename(id = sample_id, age = age_years,
         whr = whr_cm_cm, waist_cm = wc_cm,
         tnfa_ng_l = tn_fa_ng_l, adiponectin_ug_ml = adiponectin_mg_l,
         insulin_uui_ml = fasting_insulin_m_u_l) %>%
  mutate(sex = "Female",
         hip_cm = waist_cm/whr,
         chol_mg_dl = cholesterol_mmol_l * 38.67,
         hdl_mg_dl = hdl_mmol_l * 38.67,
         ldl_mg_dl = ldl_mmol_l * 38.67,
         tri_mg_dl = triglycerides_mmol_l * 88.57,
         hba1c_perc = (hb_a1c_mmol_mol/10.929) + 2.15,
         glu_mg_dl = fasting_glucose_mmol_l * 18.018) %>%
  select(-c(cholesterol_mmol_l, hdl_mmol_l, ldl_mmol_l, 
            triglycerides_mmol_l, hb_a1c_mmol_mol, fasting_glucose_mmol_l))
# lipid conversion factors: https://www.ncbi.nlm.nih.gov/books/NBK83505/
# hba1c conversion: https://pmc.ncbi.nlm.nih.gov/articles/PMC10398719/
# glucose: https://www.ncbi.nlm.nih.gov/books/NBK348987/

colnames(karlsson) %in% colnames(ai4food)
colnames(karlsson)[!colnames(karlsson) %in% colnames(ai4food)]


## ----feng--------------------------------------------------------------------------------------------------------------------------------------------
feng <- read_csv("feng_metadata_20250820.csv",
                     na = c("", "NA", "n.a", "n")) %>%
  select(subject_id, MetObesity, gender_1_male_2_female,
         age_yrs, bmi, waist_cm, hip_cm, whr_waist_hip_ratio, 
         fasting_insulin_m_u_m_l, fasting_glucose_mg_l, homa_index, hba1c_percent, 
         crp_mg_l, tg_mg_l, ldl_mg_l, hdl_mg_l) %>%
  rename(id = subject_id, sex = gender_1_male_2_female, age = age_yrs,
         bmi_kg_m2 = bmi, whr = whr_waist_hip_ratio,
         homa = homa_index, hba1c_perc = hba1c_percent,
         glu_mg_dl = fasting_glucose_mg_l,
         insulin_uui_ml = fasting_insulin_m_u_m_l,
         crp_mg_dl = crp_mg_l, tri_mg_dl = tg_mg_l,
         hdl_mg_dl = hdl_mg_l, ldl_mg_dl = ldl_mg_l) %>%
  mutate(sex = if_else(sex == 1, "Male", "Female"))

all(colnames(feng) %in% colnames(ai4food))


## ----bindRows----------------------------------------------------------------------------------------------------------------------------------------
data_export <- bind_rows(
  ai4food %>% mutate(study_name = "ai4food"),
  metacardis %>% mutate(study_name = "metacardis"), 
  karlsson %>% mutate(study_name = "karlsson"), 
  feng %>% mutate(study_name = "feng"))  %>%
  ## OUTLIERS
  mutate(hip_cm = if_else(hip_cm < 50, NA, hip_cm),
         whr = if_else(whr > 5, NA, whr)
         ) %>%
  filter(!is.na(MetObesity))


## ----physeq------------------------------------------------------------------------------------------------------------------------------------------
physeq <- readRDS("physe_adj.rds")
data_export %>% group_by(study_name) %>% count()

to.remove <- data_export %>% 
  filter(study_name == "feng", 
         !id %in% sample_names(physeq)) %>%
  pull(id)
data_export <- data_export %>% filter(!id %in% to.remove)


## ----models------------------------------------------------------------------------------------------------------------------------------------------

data$MetObesity <- factor(data$MetObesity, 
                          levels = c("MHNO", "MHO", "MUNO", "MUO"))
data$sex        <- factor(data$sex) 

var.names <- c( "hip_cm", "waist_cm", "whr", 
                "dyastolic_blood_pressure_mmhg", "systolic_blood_pressure_mmhg",
                "glu_mg_dl", "hba1c_perc", "homa", "insulin_uui_ml",
                "tri_mg_dl", "chol_mg_dl", "hdl_mg_dl", "ldl_mg_dl",
                "tnf_a_ui_ml", "crp_mg_dl", "adiponectin_ug_ml")

results <- list()

for (v in var.names) {
  
  # --- 1) Stat summary per group ---
  groups_summary <- data %>%
    group_by(MetObesity) %>%
    summarise(
      n   = sum(!is.na(.data[[v]])),
      med = median(.data[[v]], na.rm = TRUE),
      q1  = quantile(.data[[v]], 0.25, na.rm = TRUE),
      q3  = quantile(.data[[v]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(descriptive = ifelse(
      n > 0,
      paste0(round(med, 2), " (", round(q1, 2), ", ", round(q3, 2), ") [", n, "]"),
      NA_character_
    )) %>%
    select(MetObesity, descriptive)
  
  # --- 2) Adjusted linear models ---
  tmp <- data[, c(v, "MetObesity", "age", "sex", "bmi_kg_m2")]
  tmp <- tmp[complete.cases(tmp), , drop = FALSE]
  tmp$MetObesity <- droplevels(tmp$MetObesity)
  tmp$sex <- droplevels(tmp$sex)
  
  p_val <- NA_real_
  adj_r2 <- NA_real_
  emm_summary <- data.frame(MetObesity = c("MHNO", "MHO", "MUNO", "MUO"),
                            Adjusted = NA_character_)
  
  if (nrow(tmp) > 0 && nlevels(tmp$MetObesity) >= 2 && nlevels(tmp$sex) >= 2) {
    form <- as.formula(paste(v, "~ MetObesity + age + sex + bmi_kg_m2"))
    fit  <- lm(form, data = tmp)
    
    p_tab <- drop1(fit, scope = ~ MetObesity + age + sex + bmi_kg_m2, test = "F")
    p_val <- suppressWarnings(p_tab["MetObesity", "Pr(>F)"])

    adj_r2 <- summary(fit)$adj.r.squared
    
  } else {
    message(paste0("Variable '", v, "' omitted due to lack of levels."))
  }
  
  # --- 3) Get dataframe result row ---
  fila <- data.frame(
    Variable = v,
    N        = sum(!is.na(data[[v]])),
    MHNO     = groups_summary$descriptive[groups_summary$MetObesity == "MHNO"],
    MHO      = groups_summary$descriptive[groups_summary$MetObesity == "MHO"],
    MUNO     = groups_summary$descriptive[groups_summary$MetObesity == "MUNO"],
    MUO      = groups_summary$descriptive[groups_summary$MetObesity == "MUO"],

    p_value  = ifelse(is.na(p_val), NA, p_val),
    Adj_R2   = ifelse(is.na(adj_r2), NA, round(adj_r2, 3)),
    check.names = FALSE
  ) 
  
  results[[v]] <- fila
}

results_table <- bind_rows(results)
results_table$p.adj <- p.adjust(results_table$p_value, "BH")

results_table


## ----exportTable-------------------------------------------------------------------------------------------------------------------------------------
results_table_export <- results_table %>%
  mutate(p.adj = if_else(p.adj < 0.001, "< 0.001",
                         as.character(round(p.adj, 3)))) %>%
  select(-p_value) %>%
  relocate(p.adj, .before = "Adj_R2")

write_csv(results_table_export,
          "./Table_LinearModels.csv")


## ----wilcox------------------------------------------------------------------------------------------------------------------------------------------
data %>% 
  select(MetObesity, age, sex, bmi_kg_m2) %>%
  tbl_summary(by = "MetObesity") %>% 
  add_p() %>% 
  add_q(method = "BH")

data %>% kruskal_effsize(age ~ MetObesity)
data %>% kruskal_effsize(bmi_kg_m2 ~ MetObesity)
cramer_v(data$sex, data$MetObesity)


## ----emmeans-----------------------------------------------------------------------------------------------------------------------------------------
sig.vars <- results_table %>% filter(p.adj < 0.05) %>% pull(Variable)

results_posthoc <- list()

for (v in sig.vars) {
  
  tmp <- data[, c(v, "MetObesity", "age", "sex", "bmi_kg_m2")]
  tmp <- tmp[complete.cases(tmp), , drop = FALSE]
  tmp$MetObesity <- droplevels(tmp$MetObesity)
  tmp$sex <- droplevels(tmp$sex)
  
  if (nrow(tmp) > 0 && nlevels(tmp$MetObesity) >= 2 && nlevels(tmp$sex) >= 2) {
    form <- as.formula(paste(v, "~ MetObesity + age + sex + bmi_kg_m2"))
    fit  <- lm(form, data = tmp)
    
    comp <- pairs(emmeans(fit, ~ MetObesity), adjust = "tukey")
    comp_df <- as.data.frame(comp)
    
    comp_df$Variable <- v
    results_posthoc[[v]] <- comp_df
  }
}

posthoc_table <- bind_rows(results_posthoc)

posthoc_table <- posthoc_table %>%
  select(Variable, contrast, estimate, SE, df, t.ratio, p.value) %>%
  rename(
    Comparison = contrast,
    Estimate   = estimate,
    Std_Error  = SE,
    DF         = df,
    t_value    = t.ratio,
    P_value    = p.value
  ) %>%
  mutate(Feature = case_when(
    Variable == "waist_cm" ~ "Waist (cm)",
    Variable == "whr" ~ "Waist-hip ratio",
    Variable == "dyastolic_blood_pressure_mmhg" ~ "Diastolic blood pressure (mmHg)",
    Variable == "systolic_blood_pressure_mmhg" ~ "Systolic blood pressure (mmHg)",
    Variable == "glu_mg_dl" ~ "Glucose (mg/dL)",
    Variable == "hba1c_perc" ~ "HbA1c (%)",
    Variable == "homa" ~ "HOMA-IR",
    Variable == "insulin_uui_ml" ~ "Insulin (uU/mL)",
    Variable == "tri_mg_dl" ~ "Triglycerides (mg/dL)",
    Variable == "chol_mg_dl" ~ "Total cholesterol (mg/dL)",
    Variable == "hdl_mg_dl" ~ "HDL cholesterol (mg/dL)",
    Variable == "ldl_mg_dl" ~ "LDL cholesterol (mg/dL)",
  ))

posthoc_table


## ----exportPosthoc-----------------------------------------------------------------------------------------------------------------------------------
write_csv(posthoc_table %>%
            relocate(Feature, .before = "Comparison") %>%
            select(-Variable),
  "PostHocTable.csv")


## ----postHocBoxplots---------------------------------------------------------------------------------------------------------------------------------
pal <- c("#264653", "#2A9D8F", "#703d57", "#edafb8", "#e6af2e")

var.labs <- c("Waist (cm)", "Waist-hip ratio",
              "Diastolic BP (mmHg)", 
              "Systolic BP (mmHg)",
              "Glucose (mg/dL)", 
              "Hb1Ac (%)", "HOMA-IR", 
              "Insulin ( \U00B5U/mL)",
              "Triglycerides (mg/dL)", 
              "Cholesterol (mg/dL)",
              "HDL-c (mg/dL)", 
              "LDL-c (mg/dL)")
names(var.labs) <- sig.vars


stat.test <- data.frame(
  variable2 = posthoc_table$Variable,
  .y. = "value",
  group1 = gsub(" - [A-Z]*", "", posthoc_table$Comparison),
  group2 = gsub("[A-Z]* - ", "", posthoc_table$Comparison),
  p = posthoc_table$P_value
) %>%
  mutate(p.adj.signif = case_when(
    p > 0.05 ~ "ns",
    p <= 0.05 & p > 0.01 ~ "*",
    p <= 0.01 & p > 0.001 ~ "**",
    p <= 0.001 & p > 0.0001 ~ "***",
    p <= 0.0001 ~ "****"
  ))
  
ypos_list <- list()
for (v in unique(posthoc_table$Variable)) {
  
  ymax <- max(data[[v]], na.rm = TRUE)
  
  yr <- diff(range(data[[v]], na.rm = TRUE))
  
  df_v <- posthoc_table %>% filter(Variable == v)
  n <- nrow(df_v)

  base <- ymax + yr * 0.03     

  spacing <- yr * 0.08         

  df_v <- df_v %>%
    mutate(y.position = base + (0:(n-1)) * spacing)

  ypos_list[[v]] <- df_v
}



stat.test$y.position <- bind_rows(ypos_list)$y.position

plot_metadata <- data %>%
  select(MetObesity, all_of(sig.vars)) %>%
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
             labeller = labeller(variable2 = var.labs)
             ) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
 stat_pvalue_manual(
  stat.test,
  hide.ns = TRUE,
  step.increase = 0,
  tip.length = 0.01,
  bracket.shorten = 0.05,
  bracket.size = 0.3,
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
plot_metadata


## ----exportBoxplots----------------------------------------------------------------------------------------------------------------------------------
ggsave("./suppFigBoxplots.png",
       plot = plot_metadata,
       bg = "white",
       height = 13.2,
       width = 14,
       units = "cm",
       dpi = 1200)


## ----highrisk----------------------------------------------------------------------------------------------------------------------------------------
data %>% 
  filter(bmi_kg_m2 >= 40) %>% 
  group_by(MetObesity) %>% 
  count()


## ----dunnAgeBMI--------------------------------------------------------------------------------------------------------------------------------------
stat.test.2 <- data %>%
  select(MetObesity,
         age, bmi_kg_m2) %>% 
  reshape2::melt(id.vars = "MetObesity") %>% 
  mutate(variable2 = variable) %>%
  group_by(variable2) %>%
  dunn_test(value ~ MetObesity,
            p.adj = "BH") %>%
  add_xy_position(scales = "free_y")
stat.test.2


## ----exportAgeSex------------------------------------------------------------------------------------------------------------------------------------
write_csv(stat.test.2 %>%
  mutate(
    Comparison = paste(group1, "-", group2),
    Variable   = variable2,
    p.adj = if_else(p.adj < 0.001, 
                    "<0.001", 
                    as.character(round(p.adj, 2)))
  ) %>%
  select(Variable, Comparison, statistic, p, p.adj),
  "PostHocTable_AgeBMI.csv")


## ----plotAgeBMI--------------------------------------------------------------------------------------------------------------------------------------
var.labs <- c("Age", "BMI")
names(var.labs) <- c("age", "bmi_kg_m2")


plot_metadata.2 <- data %>%
  select(MetObesity, age, bmi_kg_m2) %>%
  reshape2::melt(id.vars = "MetObesity") %>%
  rename(variable2 = variable) %>%
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
  facet_wrap(. ~ variable2, 
             ncol = 3,
             scales = "free_y",
             labeller = labeller(variable2 = var.labs)) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  stat_pvalue_manual(stat.test.2,
                     hide.ns = TRUE,
                     tip.length = 0.01,
                     bracket.shorten = 0.05,
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
plot_metadata.2



## ----arrange-----------------------------------------------------------------------------------------------------------------------------------------
p <- ggarrange(plot_metadata.2,
               plot_metadata,
               nrow = 2, heights = c(1, 4))
ggsave("./suppFigBoxplots_ageBMI.png",
       plot = p,
       bg = "white",
       height = 20,
       width = 14,
       units = "cm",
       dpi = 1200)

