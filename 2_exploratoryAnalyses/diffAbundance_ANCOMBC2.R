# Analisis Abundancia Diferencial ###############
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(patchwork)

library(phyloseq)
library(microbiome)
library(ANCOMBC)
library(DT)
library(doParallel)



iniDir <- "~/ai4food"
setwd(paste0(iniDir, "/mho"))

# load physeq objects:
physeq <- readRDS("adj_physeq_MMUPHIN_mpa30_20240911.rds")
# transformar a cuentas:
physeq.abs <- transform_sample_counts(physeq, function(x) round(1E6 * x))

# ANCOM-BC W/ VIGNETTE RECOMMENDATIONS #########################################
# To control the FDR arising from multiple testing, we opt for the
# Holm-Bonferroni method over the Benjamini-Hochberg (BH) procedure,
# especially when dealing with large sample sizes where statistical power
# isn’t the primary concern. The Holm-Bonferroni method, accommodating
# any dependence structure among p-values, is known to be robust against
# inaccuracies in p-values, an issue often seen in DA analysis. Figures
# below display only results significant after the Holm-Bonferroni adjustment.

set.seed(123)
cl <- makePSOCKcluster(4)
registerDoParallel(cl)

output = ancombc2(data = physeq.abs, tax_level = "Species",
                  fix_formula = "MetObesity", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
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
                  # trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                  #                                             nrow = 2, 
                  #                                             byrow = TRUE),
                  #                                      matrix(c(-1, 0, 1, -1),
                  #                                             nrow = 2, 
                  #                                             byrow = TRUE),
                  #                                      matrix(c(1, 0, 1, -1),
                  #                                             nrow = 2, 
                  #                                             byrow = TRUE)),
                  #                      node = list(2, 2, 1),
                  #                      solver = "ECOS",
                  #                      # B = 10
                                       )

saveRDS(output, "./diffAbundance_jan25/ANCOMBC2_20250114.RDS")
stopCluster(cl)

## ANCOM-BC2 Multiple Pairwise Comparisons #####################################
# The ANCOM-BC2 methodology for multiple pairwise comparisons is designed to 
# identify taxa that exhibit differential abundance between any two groups 
# within a set of three or more experimental groups, all while maintaining 
# control over the mdFDR.
# For instance, in our analysis focusing on the categories “lean”, “overweight”, 
# and “obese”, the output provides: 1) log fold changes, 2) standard errors, 
# 3) test statistics, 4) p-values, 5) adjusted p-values, 6) indicators denoting 
# whether the taxon is differentially abundant (TRUE) or not (FALSE), and 7) 
# indicators denoting whether the taxon passed the sensitivity analysis (TRUE) 
# or not (FALSE).
# In the subsequent heatmap, each cell represents a log fold-change (in natural 
# log) value. Entries highlighted in black have successfully passed the 
# sensitivity analysis for pseudo-count addition.
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

library(RColorBrewer)
library(ggtext)
p <- df_fig_pair %>%
  mutate(taxon = gsub("_", " ", taxon)) %>%
  mutate(taxon = paste0("*", paste0(taxon, "*"))) %>%
  ggplot(aes(y = group, x = taxon, fill = value)) + 
  geom_tile(color = "gray50") +
  # colores sacados de paleta de RColorBrewer RdBu
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

ggsave("./diffAbundance_jan25/figANCOMBC2_v2.png",
       plot = p,
       bg = "white",
       height = 9.7,
       width = 15,
       units = "cm",
       dpi = 1200)
