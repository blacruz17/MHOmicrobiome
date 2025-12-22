library(tidyverse)
library(vegan)
library(gridExtra)
library(ggpubr)
library(phyloseq)
library(microbiome)
library(genefilter)
library(microViz)
library(MMUPHin)
library(ape)
library(patchwork)
library(gtsummary)
library(rstatix)
library(emmeans)

rm(list = ls())


## ----load------------------------------------------------------------------------------------------------------------
physeq <- readRDS("physeq_adj.rds")
pal <- c("#264653", "#2A9D8F", "#703d57", "#edafb8", "#e6af2e")


## ----alpha-----------------------------------------------------------------------------------------------------------
alphadiv <- alpha(physeq) %>%
  rownames_to_column(var = "id")
alphadiv_df <- sample_data(physeq) %>%
  data.frame() %>%
  rownames_to_column(var = "id") %>%
  left_join(alphadiv, by = "id")

## ----linearModel-------------------------------------------------------------
alphadiv_df$MetObesity <- factor(alphadiv_df$MetObesity, 
                          levels = c("MHNO", "MHO", "MUNO", "MUO"))
alphadiv_df$sex        <- factor(alphadiv_df$sex) 

var.names <- c("chao1", "diversity_gini_simpson", "diversity_shannon")
results <- list()

for (v in var.names) {
  
  # --- 1) Summary stats ---
  summary_groups <- alphadiv_df %>%
    group_by(MetObesity) %>%
    summarise(
      n   = sum(!is.na(.data[[v]])),
      med = median(.data[[v]], na.rm = TRUE),
      q1  = quantile(.data[[v]], 0.25, na.rm = TRUE),
      q3  = quantile(.data[[v]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(description = ifelse(
      n > 0,
      paste0(round(med, 2), " (", round(q1, 2), ", ", round(q3, 2), ") [", n, "]"),
      NA_character_
    )) %>%
    select(MetObesity, description)
  
  # --- 2) Adj linear model ---
  tmp <- alphadiv_df[, c(v, "MetObesity", "age", "sex", "bmi_kg_m2")]
  tmp <- tmp[complete.cases(tmp), , drop = FALSE]
  tmp$MetObesity <- droplevels(tmp$MetObesity)
  tmp$sex <- droplevels(tmp$sex)
  tmp$study_name <- droplevels(tmp$study_name)
  
  p_val <- NA_real_
  adj_r2 <- NA_real_
  emm_resumen <- data.frame(MetObesity = c("MHNO", "MHO", "MUNO", "MUO"),
                            Adjusted = NA_character_)
  
  if (nrow(tmp) > 0 && nlevels(tmp$MetObesity) >= 2 && nlevels(tmp$sex) >= 2) {
    form <- as.formula(paste(v, "~ MetObesity + age + sex + bmi_kg_m2"))
    fit  <- lm(form, data = tmp)
    
    # p-valor global del grupo (test parcial con drop1)
    p_tab <- drop1(fit, scope = ~ MetObesity + age + sex + bmi_kg_m2, test = "F")
    p_val <- suppressWarnings(p_tab["MetObesity", "Pr(>F)"])
    
    # R² Adjusted
    adj_r2 <- summary(fit)$adj.r.squared
    
  } else {
    message(paste0("Variable '", v, "' omitted due to lack of levels."))
  }
  
  # --- 3) Get dataframe row w/ results ---
  fila <- data.frame(
    Variable = v,
    N        = sum(!is.na(alphadiv_df[[v]])),
    MHNO     = summary_groups$description[summary_groups$MetObesity == "MHNO"],
    MHO      = summary_groups$description[summary_groups$MetObesity == "MHO"],
    MUNO     = summary_groups$description[summary_groups$MetObesity == "MUNO"],
    MUO      = summary_groups$description[summary_groups$MetObesity == "MUO"],
    p_value  = ifelse(is.na(p_val), NA, p_val),
    Adj_R2   = ifelse(is.na(adj_r2), NA, round(adj_r2, 3)),
    check.names = FALSE
  ) 
  
  results[[v]] <- fila
}

results_table <- bind_rows(results)
results_table$p.adj <- p.adjust(results_table$p_value, "BH")

results_table

results_table_export <- results_table %>%
  mutate(p.adj = if_else(p.adj < 0.001, "<0.001",
                         as.character(round(p.adj, 3)))) %>%
  select(-p_value) %>%
  relocate(p.adj, .before = "Adj_R2")
results_table_export

write_csv(results_table_export,
          "./Table_AlphaDiv_LinearModels.csv")

## ----PostHocTukey------------------------------------------------------------
sig.vars <- results_table %>% filter(p.adj < 0.05) %>% pull(Variable)

results_posthoc <- list()

for (v in sig.vars) {
  
  tmp <- alphadiv_df[, c(v, "MetObesity", "age", "sex", "bmi_kg_m2")]
  tmp <- tmp[complete.cases(tmp), , drop = FALSE]
  tmp$MetObesity <- droplevels(tmp$MetObesity)
  tmp$sex <- droplevels(tmp$sex)
  tmp$study_name <- droplevels(tmp$study_name)
  
  if (nrow(tmp) > 0 && nlevels(tmp$MetObesity) >= 2 && nlevels(tmp$sex) >= 2) {
    form <- as.formula(paste(v, "~ MetObesity + age + sex + bmi_kg_m2"))
    fit  <- lm(form, data = tmp)
    
    comp <- pairs(emmeans(fit, ~ MetObesity), adjust = "tukey")
    comp_df <- as.data.frame(comp)
  
    comp_df$Variable <- v
    results_posthoc[[v]] <- comp_df
  }
}

posthoc_table <- bind_rows(results_posthoc) %>%
  select(Variable, contrast, estimate, SE, df, t.ratio, p.value) %>%
  rename(
    Comparison = contrast,
    Estimate   = estimate,
    Std_Error  = SE,
    DF         = df,
    t_value    = t.ratio,
    P_value    = p.value
  ) %>%
  mutate(P_value = if_else(P_value < 0.001, "<0.001",
                           as.character(round(P_value, 3))),
         Feature = case_when(
           Variable == "chao1" ~ "Chao1",
           Variable == "diversity_shannon" ~ "Shannon"
         ))

posthoc_table 

write_csv(posthoc_table %>%
            relocate(Feature, .before = "Comparison") %>%
            select(-Variable),
          "PostHocTable_AlphaDiv.csv")

## ----alphaplots------------------------------------------------------------------------------------------------------
pal <- c("#264653", "#2A9D8F", "#703d57", "#edafb8", "#e6af2e")

var.labs <- c("Chao1", "Shannon", "Simpson")
names(var.labs) <- c("chao1", "diversity_shannon", "diversity_gini_simpson")
var.labs

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
  
  ymax <- max(alphadiv_df[[v]], na.rm = TRUE)
  
  df_v <- posthoc_table %>% filter(Variable == v)
  
  n <- nrow(df_v)
  margin <- (ymax * 0.05) 
  df_v <- df_v %>%
    mutate(y.position = ymax + (1:n) * margin)
  
  ypos_list[[v]] <- df_v
}

stat.test$y.position <- bind_rows(ypos_list)$y.position

alpha.p <- alphadiv_df  %>% 
  select(MetObesity, 
         chao1, 
         diversity_gini_simpson,
         diversity_shannon) %>%
  reshape2::melt(id.vars = "MetObesity") %>%
  rename(variable2 = variable) %>%
  ggplot(aes(x = MetObesity, y = value))+
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
                     name = "")  +
  facet_wrap(. ~ variable2, 
             ncol = 3,
             scales = "free_y",
             labeller = labeller(variable2 = var.labs))  + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  stat_pvalue_manual(stat.test,
                     hide.ns = TRUE,
                     tip.length = 0.01,
                     bracket.shorten = 0.05,
                     step.increase = .02,
                     bracket.size = 0.3,
                     label.size = 4) +
  theme_linedraw() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    plot.title = element_text(face = 'bold'),
    axis.text = element_text(size = 12),
    strip.text = element_text(face = 'bold', size = 14,
                              color = "black"),
    strip.background = element_rect(
      color="white", fill="white",  linetype="solid"
    ))
alpha.p



## ----saveAlphaPlot---------------------------------------------------------------------------------------------------
ggsave("./alphaDiv_Adjusted.png", 
       plot = alpha.p,
       width = 8, height = 3.15, units = "in", dpi = 1200)


## ----permanovaMetHealth----------------------------------------------------------------------------------------------
physeq.clr <- microbiome::transform(physeq, "clr")
aitchison_dist <- distance(physeq.clr, method = "euclidean")
aitchison_adonis <- adonis2(aitchison_dist ~ MetObesity,
                            data = data.frame(sample_data(physeq.clr)))

ait_adonis_adj <- adonis2(aitchison_dist ~ MetObesity  + age + sex + bmi_kg_m2, 
                          data = data.frame(sample_data(physeq.clr)), 
                          by = "terms")

ait_adonis_adj2 <- adonis2(aitchison_dist ~  age + sex + bmi_kg_m2, 
                          data = data.frame(sample_data(physeq.clr)), 
                          by = "terms")
## ----showAdonis------------------------------------------------------------------------------------------------------
ait_adonis_adj; ait_adonis_adj2


## ----plotAitchi------------------------------------------------------------------------------------------------------
aitchison_pcoa  <- ordinate(physeq.clr, "PCoA", distance = aitchison_dist)

p1 <- plot_ordination(physeq.clr,
                      aitchison_pcoa,
                      color = "MetObesity") +
  geom_point(size = .6) +
  scale_colour_manual(values = pal) +
            theme_linedraw() +
            stat_ellipse(aes(color = MetObesity),
                         linewidth = .4) +
            theme(legend.position = "none")
p1$layers <- p1$layers[-1] # keep only smaller points

p2 <- plot_ordination(physeq.clr,
                      aitchison_pcoa,
                      color = "MetObesity",
                      axes = c(1,3)) +
  geom_point(size = .6) +
  scale_colour_manual(values = pal) +
  theme_linedraw() +
  stat_ellipse(aes(color = MetObesity),
               linewidth = .4) +
  theme(legend.position = "none")
p2$layers <- p2$layers[-1]

p3 <- plot_ordination(physeq.clr,
                      aitchison_pcoa,
                      color = "MetObesity",
                      axes = c(2,3)) +
  geom_point(size = .6) +
  scale_colour_manual(values = pal) +
  theme_linedraw() +
  stat_ellipse(aes(color = MetObesity),
               linewidth = .4) +
  theme(legend.title = element_blank())
p3$layers <- p3$layers[-1]

beta.p <- p1 + p2 + p3 & theme(axis.text = element_text(size = 10),
                               axis.title = element_text(size = 12),
                               legend.text = element_text(size = 10),
                               panel.border = element_rect(colour = "black",
                                                           fill=NA, 
                                                           linewidth=.3))


## ----saveBetaPlot----------------------------------------------------------------------------------------------------
beta.p
ggsave("./betaDiv.png",
       width = 8, height = 2.65, units = "in", dpi = 1200)
