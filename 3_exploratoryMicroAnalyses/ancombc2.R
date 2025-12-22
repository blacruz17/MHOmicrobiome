library(tidyverse)
library(gridExtra)
library(ggpubr)
library(patchwork)
library(phyloseq)
library(microbiome)
library(ANCOMBC)
library(DT)
library(doParallel)
library(RColorBrewer)
library(ggtext)

rm(list = ls())

## ----load------------------------------------------------------------------------------------------------------------
physeq <- readRDS("physeq_adj.rds")
pal <- c("#264653", "#2A9D8F", "#703d57", "#edafb8", "#e6af2e")
physeq.abs <- transform_sample_counts(physeq, function(x) round(1E6 * x))

set.seed(123)
cl <- makePSOCKcluster(4)
registerDoParallel(cl)

output = ancombc2(data = physeq.abs, 
                  tax_level = "Species",
                  fix_formula = "MetObesity + age + sex + bmi_kg_m2", 
                  rand_formula = NULL,
                  p_adj_method = "holm", 
                  pseudo_sens = TRUE,
                  prv_cut = 0, 
                  lib_cut = 1000, s0_perc = 0.05,
                  group = "MetObesity", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = TRUE, pairwise = TRUE, dunnet = TRUE, 
                  trend = FALSE,
                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
)

stopCluster(cl)

res_pair = output$res_pair

df_fig_pair1 = res_pair %>%
  dplyr::filter(diff_MetObesityMHO == 1 |
                  diff_MetObesityMUO == 1 | 
                  diff_MetObesityMUNO == 1 |
                  diff_MetObesityMUNO_MetObesityMHO == 1 |
                  diff_MetObesityMUO_MetObesityMHO == 1 |
                  diff_MetObesityMUO_MetObesityMUNO == 1
  ) %>%
  dplyr::mutate(lfc1 = ifelse(diff_MetObesityMHO == 1, 
                              round(lfc_MetObesityMHO, 2), 0),
                lfc2 = ifelse(diff_MetObesityMUO == 1, 
                              round(lfc_MetObesityMUO, 2), 0),
                lfc3 = ifelse(diff_MetObesityMUNO == 1, 
                              round(lfc_MetObesityMUNO, 2), 0),
                lfc4 = ifelse(diff_MetObesityMUNO_MetObesityMHO == 1, 
                              round(lfc_MetObesityMUNO_MetObesityMHO, 2), 0),
                lfc5 = ifelse(diff_MetObesityMUO_MetObesityMHO == 1, 
                              round(lfc_MetObesityMUO_MetObesityMHO, 2), 0),
                lfc6 = ifelse(diff_MetObesityMUO_MetObesityMUNO == 1, 
                              round(lfc_MetObesityMUO_MetObesityMUNO, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc6, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_pair2 = res_pair %>%
  dplyr::filter(diff_MetObesityMHO == 1 |
                  diff_MetObesityMUO == 1 | 
                  diff_MetObesityMUNO == 1 |
                  diff_MetObesityMUNO_MetObesityMHO == 1 |
                  diff_MetObesityMUO_MetObesityMHO == 1 |
                  diff_MetObesityMUO_MetObesityMUNO == 1
  ) %>%
  dplyr::mutate(lfc1 = "black",
                lfc2 = "black",
                lfc3 = "black",
                lfc4 = "black",
                lfc5 = "black",
                lfc6 = "black") %>%
  tidyr::pivot_longer(cols = lfc1:lfc6, 
                      names_to = "group", values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_pair = df_fig_pair1 %>%
  dplyr::left_join(df_fig_pair2, by = c("taxon", "group"))

df_fig_pair$group = recode(df_fig_pair$group, 
                           `lfc1` = "MHO - MHNO",
                           `lfc2` = "MUO - MHNO",
                           `lfc3` = "MUNO - MHNO",
                           `lfc4` = "MUNO - MHO",
                           `lfc5` = "MUO - MHO",
                           `lfc6` = "MUO - MUNO")
df_fig_pair$group = factor(df_fig_pair$group,
                           levels = c("MHO - MHNO",
                                      "MUO - MHO",
                                      "MUNO - MHO",
                                      "MUO - MHNO",
                                      "MUNO - MHNO",
                                      "MUO - MUNO"))

lo = floor(min(df_fig_pair$value))
up = ceiling(max(df_fig_pair$value))
mid = (lo + up)/2
mid = 0

p <- df_fig_pair %>%
  mutate(taxon = gsub("_", " ", taxon)) %>%
  mutate(taxon = paste0("*", paste0(taxon, "*"))) %>%
  ggplot(aes(y = group, x = taxon, fill = value)) + 
  geom_tile(color = "gray50") +
  scale_fill_gradient2(low = "#313695", high = "#A50026" , mid = "white" , 
                       na.value = "white" , midpoint = mid, limit = c(lo, up),
                       name = "LFC") +
  geom_text(aes(taxon, group, label = value, color = color), size = 2) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_bw() +
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1),
        axis.text = element_text(color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))

p

ggsave("./diffAbundance_adjusted.png",
       plot = p,
       bg = "white",
       height = 9.7,
       width = 15,
       units = "cm",
       dpi = 1200)
