#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Plotting utilities for ANCOM-BC2 pairwise results (GO terms)
#
# This script loads an ANCOM-BC2 result object and produces:
# (1) A heatmap of log-fold changes (LFC) across pairwise MetObesity contrasts.
# (2) A facetted forest plot for a curated set of GO terms across all contrasts.
# (3) Two single-contrast forest plots (MUNO vs MHNO, MUO vs MHNO) with p-values.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(tidyverse)   
  library(grid)
  library(gridExtra)
})

# Load ANCOM-BC2 results

out <- readRDS("Ancom-bc/diffAbundance_sep25_GO/ANCOMBC2_20251013_allvariable_top1000.RDS")
res_pair <- out$res_pair

# (1) Heatmap of LFCs for all significant pairwise contrasts

df_fig_pair1 <- res_pair %>%
  filter(
    diff_MetObesityMHO == 1 |
      diff_MetObesityMUO == 1 |
      diff_MetObesityMUNO == 1 |
      diff_MetObesityMUNO_MetObesityMHO == 1 |
      diff_MetObesityMUO_MetObesityMHO == 1 |
      diff_MetObesityMUO_MetObesityMUNO == 1
  ) %>%
  mutate(
    lfc1 = ifelse(diff_MetObesityMHO == 1, round(lfc_MetObesityMHO, 2), 0),
    lfc2 = ifelse(diff_MetObesityMUO == 1, round(lfc_MetObesityMUO, 2), 0),
    lfc3 = ifelse(diff_MetObesityMUNO == 1, round(lfc_MetObesityMUNO, 2), 0),
    lfc4 = ifelse(diff_MetObesityMUNO_MetObesityMHO == 1, round(lfc_MetObesityMUNO_MetObesityMHO, 2), 0),
    lfc5 = ifelse(diff_MetObesityMUO_MetObesityMHO == 1, round(lfc_MetObesityMUO_MetObesityMHO, 2), 0),
    lfc6 = ifelse(diff_MetObesityMUO_MetObesityMUNO == 1, round(lfc_MetObesityMUO_MetObesityMUNO, 2), 0)
  ) %>%
  pivot_longer(cols = lfc1:lfc6, names_to = "group", values_to = "value") %>%
  arrange(taxon)

df_fig_pair2 <- res_pair %>%
  filter(
    diff_MetObesityMHO == 1 |
      diff_MetObesityMUO == 1 |
      diff_MetObesityMUNO == 1 |
      diff_MetObesityMUNO_MetObesityMHO == 1 |
      diff_MetObesityMUO_MetObesityMHO == 1 |
      diff_MetObesityMUO_MetObesityMUNO == 1
  ) %>%
  mutate(
    lfc1 = ifelse(passed_ss_MetObesityMHO == 1 & diff_MetObesityMHO == 1, "black", "gray50"),
    lfc2 = ifelse(passed_ss_MetObesityMUO == 1 & diff_MetObesityMUO == 1, "black", "gray50"),
    lfc3 = ifelse(passed_ss_MetObesityMUNO == 1 & diff_MetObesityMUNO == 1, "black", "gray50"),
    lfc4 = ifelse(passed_ss_MetObesityMUNO_MetObesityMHO == 1 & diff_MetObesityMUNO_MetObesityMHO == 1, "black", "gray50"),
    lfc5 = ifelse(passed_ss_MetObesityMUO_MetObesityMHO == 1 & diff_MetObesityMUO_MetObesityMHO == 1, "black", "gray50"),
    lfc6 = ifelse(passed_ss_MetObesityMUO_MetObesityMUNO == 1 & diff_MetObesityMUO_MetObesityMUNO == 1, "black", "gray50")
  ) %>%
  pivot_longer(cols = lfc1:lfc6, names_to = "group", values_to = "color") %>%
  arrange(taxon)

df_fig_pair <- df_fig_pair1 %>%
  left_join(df_fig_pair2, by = c("taxon", "group"))

# Provide human-readable labels and an explicit ordering of contrasts.
df_fig_pair$group <- recode(
  df_fig_pair$group,
  `lfc1` = "MHO - MHNO",
  `lfc2` = "MUO - MHNO",
  `lfc3` = "MUNO - MHNO",
  `lfc4` = "MUNO - MHO",
  `lfc5` = "MUO - MHO",
  `lfc6` = "MUO - MUNO"
)
df_fig_pair$group <- factor(
  df_fig_pair$group,
  levels = c("MHO - MHNO", "MUO - MHO", "MUNO - MHO", "MUO - MHNO", "MUNO - MHNO", "MUO - MUNO")
)

# Use symmetric midpoint at 0 to emphasise directionality in LFC.
lo <- floor(min(df_fig_pair$value))
up <- ceiling(max(df_fig_pair$value))
mid <- 0

fig_pair <- df_fig_pair %>%
  ggplot(aes(x = group, y = taxon, fill = value)) +
  geom_tile(color = "gray50") +
  scale_fill_gradient2(
    low = "#313695", high = "#A50026", mid = "white",
    na.value = "white", midpoint = mid, limit = c(lo, up), name = NULL
  ) +
  geom_text(aes(group, taxon, label = value, color = color), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold changes in pairwise comparisons") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5))

print(fig_pair)

ggsave("Ancom-bc/diffAbundance_sep25_GO/fig_pair_allvariable_top1000.pdf",plot = fig_pair, width = 8, height = 6,, dpi = 800)

# (2) Facetted forest plot for a curated set of GO terms across all contrasts

go_interes <- c("GO:0004143", "GO:0004574", "GO:0004617", "GO:0050112", "GO:0004309", "GO:0052692", "GO:0070626")

go_names <- c(
  "GO:0004143" = "TP-dependent diacylglycerol kinase activity",
  "GO:0004574" = "oligo-1,6-glucosidase activity",
  "GO:0004617" = "phosphoglycerate dehydrogenase activity",
  "GO:0050112" = "inositol 2-dehydrogenase (NAD+) activity",
  "GO:0004309" = "exopolyphosphatase activity",
  "GO:0052692" = "raffinose alpha-galactosidase activity",
  "GO:0070626" = "(S)-2-(5-amino-1-(5-phospho-D-ribosyl)imidazole-4-carboxamido) succinate lyase (fumarate-forming) activity"
)

df_go <- res_pair %>% filter(taxon %in% go_interes)

# Reshape to long format, aligning LFC, SE, q-values, and robust-differential flags.
df_long <- df_go %>%
  select(taxon,
         starts_with("lfc_"),
         starts_with("se_"),
         starts_with("p_"),       
         starts_with("q_"),
         starts_with("diff_robust_")) %>%
  pivot_longer(
    cols = -taxon,
    names_to = c(".value", "comparison"),
    #names_pattern = "(lfc|se|q|diff_robust)_(.*)"
    names_pattern = "(lfc|se|p|q|diff_robust)_(.*)"  
  ) %>%
  mutate(
    comparison = recode(comparison,
      "MetObesityMHO" = "MHO - MHNO",
      "MetObesityMUO" = "MUO - MHNO",
      "MetObesityMUNO" = "MUNO - MHNO",
      "MetObesityMUNO_MetObesityMHO" = "MUNO - MHO",
      "MetObesityMUO_MetObesityMHO"  = "MUO - MHO",
      "MetObesityMUO_MetObesityMUNO" = "MUO - MUNO"
    ),
    #robust = diff_robust == 1,
    CI_low = lfc - 1.96 * se,
    CI_high = lfc + 1.96 * se,
    #color = ifelse(robust, "black", "gray60")
    sig = q < 0.05,
    #sig = p < 0.05,
    color = ifelse(sig, "black", "gray60"),
     # etiqueta de q solo si es significativo y es MUO-MHNO o MUNO-MHNO
    label_q = ifelse(
      sig & comparison %in% c("MUO - MHNO", "MUNO - MHNO"),
      paste0("q=", signif(q, 2)),   # o 3 si quieres más decimales
      NA_character_
    )

  )
comparison_levels <- c("MHO - MHNO", "MUO - MHNO", "MUNO - MHNO", "MUNO - MHO", "MUO - MHO", "MUO - MUNO")
go_levels <- rev(go_interes)

p_go <- ggplot(df_long, aes(x = lfc, y = factor(taxon, levels = go_levels))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, color = color),
                 height = 0.15, linewidth = 0.6) +
  geom_point(aes(color = color), size = 2) +
  geom_text(
    aes(x = CI_high + 0.03,      
        label = label_q),
    hjust = 0,                  
    size  = 3,
    na.rm = TRUE
  ) +
  facet_wrap(~comparison, ncol = 2) +
  scale_color_identity() +
  theme_bw(base_size = 12) +
  theme(
    strip.background   = element_rect(fill = "grey90"),
    axis.title.y       = element_blank(),
    axis.title.x       = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey95"),
    axis.text.y        = element_text(size = 10)
  ) +
  labs(
    x = "Log fold change (LFC)"
  ) +
  coord_cartesian(xlim = c(min(df_long$CI_low),
                           max(df_long$CI_high) + 0.15))

ggsave("ANCOMBC2_forest.png",  plot = p_go, width = 8, height = 5, dpi = 800)

# (3) Two single-contrast forest plots with p-values (MUNO vs MHNO; MUO vs MHNO)

# Keep only comparisons against the MHNO reference.
comparaciones_keep <- c("MetObesityMUNO", "MetObesityMUO")

df_long2 <- df_go %>%
  select(
    taxon,
    any_of(paste0("lfc_", comparaciones_keep)),
    any_of(paste0("se_",  comparaciones_keep)),
    any_of(paste0("p_",   comparaciones_keep)),
    any_of(paste0("diff_robust_", comparaciones_keep))
  ) %>%
  pivot_longer(
    cols = -taxon,
    names_to = c(".value", "comparison_raw"),
    names_pattern = "(lfc|se|p|diff_robust)_(.*)"
  ) %>%
  filter(comparison_raw %in% comparaciones_keep) %>%
  mutate(
    comparison = recode(
      comparison_raw,
      "MetObesityMUNO" = "MUNO vs MHNO",
      "MetObesityMUO"  = "MUO vs MHNO"
    ),
    robust    = diff_robust == 1,
    CI_low    = lfc - 1.96 * se,
    CI_high   = lfc + 1.96 * se,
    color_txt = ifelse(robust, "black", "gray60")
  )

go_levels <- rev(go_interes)

plot_forest <- function(df_comp, comp_label, accent_hex) {
  dfc <- df_comp %>% filter(comparison == comp_label)
  n_rows <- length(unique(dfc$taxon))
  
  ggplot(dfc, aes(x = lfc, y = factor(taxon, levels = go_levels))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.18, linewidth = 0.6) +
    geom_point(size = 2.2) +
    geom_text(
      aes(
        label = paste0("p=", formatC(p, format = "g", digits = 3)),
        color = color_txt
      ),
      x = Inf, hjust = -0.1, size = 3.3, show.legend = FALSE, na.rm = TRUE
    ) +
    annotate("text", x = Inf, y = n_rows + 0.7, label = "p-value", hjust = -0.1, fontface = "bold") +
    scale_y_discrete(labels = identity) +
    scale_color_identity() +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.30))) +
    coord_cartesian(clip = "off") +
    labs(
      x = "Log fold change (LFC)", y = NULL,
      title = paste0("ANCOM-BC2 — ", comp_label)
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", color = accent_hex),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_line(color = "grey95"),
      plot.margin = margin(10, 90, 10, 10),
      axis.title.x = element_text(face = "bold")
    )
}

col_MUNO <- "#703d57"
col_MUO  <- "#edafb8"

p_muno <- plot_forest(df_long2, "MUNO vs MHNO", accent_hex = col_MUNO)
p_muo  <- plot_forest(df_long2, "MUO vs MHNO",  accent_hex = col_MUO)

grid.arrange(p_muno, p_muo, ncol = 1)

ggsave("ANCOMBC2_forest_MUO_vs_MHNO.png",  p_muo,  width = 8, height = 5, dpi = 800)
