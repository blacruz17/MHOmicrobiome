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

## ----load------------------------------------------------------------------------------------------------------------
physeq <- readRDS("..data/physeqMHO.rds")
pal <- c("#264653", "#2A9D8F", "#703d57", "#edafb8", "#e6af2e")


## ----alpha-----------------------------------------------------------------------------------------------------------
alphadiv <- alpha(physeq) %>%
  rownames_to_column(var = "id")
alphadiv_df <- sample_data(physeq) %>%
  data.frame() %>%
  rownames_to_column(var = "id") %>%
  left_join(alphadiv, by = "id")


## ----alphaStats------------------------------------------------------------------------------------------------------
alphadiv_df %>% 
  select(MetObesity, 
         chao1, diversity_gini_simpson,
         diversity_shannon,
         dominance_absolute,
         rarity_low_abundance) %>% 
  tbl_summary(by = MetObesity) %>% 
  add_p() %>% 
  add_q()


## ----alphaDunn-------------------------------------------------------------------------------------------------------
stat.test <- alphadiv_df %>% 
  select(MetObesity, 
         chao1, diversity_gini_simpson,
         diversity_shannon) %>% 
  reshape2::melt(id.vars = "MetObesity") %>% 
  mutate(variable2 = variable) %>%
  group_by(variable2) %>%
  dunn_test(value ~ MetObesity,
            p.adj = "BH") %>%
  add_xy_position(scales = "free_y")

stat.test


## ----alphaplots------------------------------------------------------------------------------------------------------
var.labs <- c("Chao1", "Simpson", "Shannon")
names(var.labs) <- c("chao1", 
                     "diversity_gini_simpson",
                     "diversity_shannon")
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
alpha.p
ggsave("../figures/alphaDiv.png",
       width = 8, height = 3.15, units = "in", dpi = 1200)


## ----permanovaMetHealth----------------------------------------------------------------------------------------------
physeq.clr <- microbiome::transform(physeq, "clr")
aitchison_dist <- distance(physeq.clr, method = "euclidean")
aitchison_adonis <- adonis2(aitchison_dist ~ MetObesity,
                                data = data.frame(sample_data(physeq.clr)))


## ----showAdonis------------------------------------------------------------------------------------------------------
aitchison_adonis


## ----permanovaStudy--------------------------------------------------------------------------------------------------
aitchison_adonis_study <- adonis2(aitchison_dist ~ study_name,
                                  data = data.frame(sample_data(physeq.clr)))
aitchison_adonis_study


## ----otherDistances--------------------------------------------------------------------------------------------------
ps.unifrac <- physeq %>%
  dist_calc("unifrac") 
ps.wunifrac <- physeq %>%
  dist_calc("wunifrac") 
ps.jaccard <- physeq %>%
  dist_calc("jaccard") 
ps.bray <- physeq %>%
  dist_calc("bray")


## ----otherPermanovas-------------------------------------------------------------------------------------------------
PERM.UF <- ps.unifrac %>% 
  dist_permanova(variables = c("MetObesity", "study_name"),
                 n_processes = 4)

PERM.WUF <- ps.wunifrac %>% 
  dist_permanova(variables = c("MetObesity", "study_name"),
                 n_processes = 4)

PERM.B <- ps.bray %>% 
  dist_permanova(variables = c("MetObesity", "study_name"),
                 n_processes = 4)

PERM.J <- ps.jaccard %>% 
  dist_permanova(variables = c("MetObesity", "study_name"),
                 n_processes = 4)


## ----otherPermanovas_show--------------------------------------------------------------------------------------------
PERM.UF
PERM.WUF
PERM.B
PERM.J


## ----plotAitchi------------------------------------------------------------------------------------------------------
p1 <- physeq %>%
  tax_transform(trans = "clr", rank = "Species") %>% 
  dist_calc("euclidean") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "MetObesity",
           # alpha = .9,
           auto_caption = NA,
           axes = 1:2,
           size = .6
           ) +
  scale_colour_manual(values = pal) +
  theme_linedraw() +
  stat_ellipse(aes(color = MetObesity),
               linetype = "dashed", linewidth = .8) +
  labs(x = "Dim1 (8.0%)",
       y = "Dim2 (4.6%)") +
  theme(legend.position = "none")

p2 <- physeq %>%
  tax_transform(trans = "clr", rank = "Species") %>% 
  dist_calc("euclidean") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "MetObesity",
           # alpha = .9,
          size = .6,
           auto_caption = NA,
           axes = c(1,3)
  ) +
  scale_colour_manual(values = pal) +
  theme_linedraw() +
  stat_ellipse(aes(color = MetObesity),
               linetype = "dashed", linewidth = .8) +
  labs(x = "Dim1 (8.0%)",
       y = "Dim3 (3.9%)") +
  theme(legend.position = "none")

p3 <- physeq %>%
  tax_transform(trans = "clr", rank = "Species") %>% 
  dist_calc("euclidean") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "MetObesity",
           # alpha = .9,
          size = .6,
           auto_caption = NA,
           axes = c(2,3)
  ) +
  scale_colour_manual(values = pal) +
  theme_linedraw() +
  stat_ellipse(aes(color = MetObesity),
               linetype = "dashed", linewidth = .8) +
  labs(x ="Dim2 (4.6%)",
       y = "Dim3 (3.9%)") +
  theme(legend.title = element_blank())

beta.p <- p1 + p2 + p3 & theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14),
                     legend.text = element_text(size = 12))


## ----saveBetaPlot----------------------------------------------------------------------------------------------------
beta.p
ggsave("../figures/betaDiv.png",
       width = 8, height = 2.65, units = "in", dpi = 1200)