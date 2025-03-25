rm(list = ls())

library(phyloseq)
library(tidyverse)
library(janitor)
library(readxl)
library(ComplexHeatmap)
library(gtsummary)
library(RColorBrewer)
library(rstatix)
library(ggpubr)
library(gt)


physeq <- readRDS("../data/physeqMHO.rds")  
data <- data.frame(sample_data(physeq))


tbl_all <- data %>%
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
tbl_all 

gt::gtsave(as_gt(tbl_all), filename = "../figures/table1.html")

dunn.tb <- data %>%
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


stat.test <- data %>%
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
  add_xy_position(scales = "free_y") %>%
  mutate(p.adj.signif = if_else(p.adj.signif == "****", 
                                "***", p.adj.signif)) # only below 0.001


pal <- c("#264653", "#2A9D8F", "#703d57", "#edafb8", "#e6af2e")

var.labs <- c("Age", "BMI", 
              "Weight", "Hip", "Waist", "Waist-hip ratio",
              "Glucose", "Hb1Ac", "HOMA-IR", "Insulin",
              "Triglycerides", "HDL", 
              "Diastolic blood pressure", "Systolic blood pressure",
              "TNF- \U03B1 (U/mL)", 
              "Adiponectin", "C-reactive protein")
names(var.labs) <- c("age", "bmi", 
         "weight", "hip", "waist", "whr",
         "glucose", "hba1c", "homa", "insulin",
         "triglycerides", "hdl", 
         "diastolic", "systolic",
         "tnfa_1",
         "adiponectin", "crp")

p.metadata_exploration <- data %>%
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
                     label.size = 1.5
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

ggsave("./figures/suppFigBoxplots.png",
       plot = p.metadata_exploration,
       bg = "white",
       height = 20,
       width = 14,
       units = "cm",
       dpi = 1200)
