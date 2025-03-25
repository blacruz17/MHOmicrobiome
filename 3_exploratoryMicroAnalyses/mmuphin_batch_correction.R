# CRAN
library(tidyverse)
library(vegan)
library(gridExtra)
library(ggpubr)

# Bioconductor
library(phyloseq)
library(microbiome)
library(genefilter)
library(microViz)
library(MMUPHin)
library(ape)

## ----import-----------------------------------------------------------------------------------------------
rm(list = ls())
physeq <- readRDS("../data/physeqMHO_raw.rds")
metaphlanToPhyloseq <- source('./metaphlanToPhyloseq.R')$value


## ----addMetadata------------------------------------------------------------------------------------------
sample_data(physeq) <- sample_data(physeq) %>% 
  data.frame() %>%
  mutate(isCurated = if_else(study_name == "AI4Food", "No", "Yes"))


## ----filter-----------------------------------------------------------------------------------------------
flist_1 <- filterfun(kOverA(k = 10, A = 0.01))
(physeq.fil <- filter_taxa(physeq = physeq,
                           flist = flist_1,
                           prune = TRUE)) 

# remove sample with all abundances equal to 0
physeq.fil <- prune_samples(sample_names(physeq.fil) != 
                              "HMP2_J05379_M_ST_T0_B0_0120_ZOZOW1T-44_H15WVBG", 
                            physeq.fil)


## ----postFilter-------------------------------------------------------------------------------------------
summary(colSums(otu_table(physeq.fil))) ; hist(colSums(otu_table(physeq.fil)))


## ----aitchison--------------------------------------------------------------------------------------------
# aitchison:
physeq.clr <- microbiome::transform(physeq.fil, "clr")

aitchison_dist <- distance(physeq.clr, method = "euclidean")
aitchison_pcoa  <- ordinate(physeq.clr, "PCoA", distance = aitchison_dist)

clr.plot <- plot_ordination(physeq.clr,
                            aitchison_pcoa, 
                            color = "study_name") + 
  geom_point() +
  theme_bw()

clr.plot


## ----microViz---------------------------------------------------------------------------------------------
## i - PCA
physeq.fil %>%
  # tax_transform("clr", rank = "Species") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>%
  ord_plot(color = "study_name", 
           shape = "MetObesity", 
           size = 2) +
  scale_colour_brewer(palette = "Set2")+
  scale_colour_brewer(palette = "Set2", aesthetics = c("fill", "colour"), name = "Study") +
  theme_bw() +
  ggside::geom_xsidedensity(aes(fill = study_name), 
                            alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = study_name), 
                            alpha = 0.5, show.legend = FALSE) +
  ggside::theme_ggside_void()

## ii - PCoA
physeq.fil %>%
  tax_transform(trans = "clr", rank = "Species") %>% # don't transform!
  dist_calc("euclidean") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "study_name", shape = "MetObesity", size = 2) +
  scale_colour_brewer(palette = "Set2") +
  scale_colour_brewer(palette = "Set2", aesthetics = c("fill", "colour"), name = "Study") +
  theme_bw() +
  ggside::geom_xsidedensity(aes(fill = study_name), 
                            alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = study_name), 
                            alpha = 0.5, show.legend = FALSE) +
  ggside::theme_ggside_void()


## ----AdonisStudyName--------------------------------------------------------------------------------------
aitchison_adonis_study <- vegan::adonis(aitchison_dist ~ study_name,
                                        data = data.frame(sample_data(physeq.clr)))
aitchison_adonis_study$aov.tab 


## ----AdonisCurated----------------------------------------------------------------------------------------
aitchison_adonis_curated <- vegan::adonis(aitchison_dist ~ isCurated,
                                        data = data.frame(sample_data(physeq.clr)))
aitchison_adonis_curated$aov.tab 


## ----mmuphin----------------------------------------------------------------------------------------------
# converts to numbers 0-1:
merged_otu <- otu_table(physeq.fil) %>% 
  apply(2, function(x) x/100)
merged_metadata <- data.frame(sample_data(physeq.fil))

# batch effect adjustment:
fit_adjust_batch <- adjust_batch(feature_abd = merged_otu,
                                 batch = "study_name",
                                 covariates = "MetObesity",
                                 data = merged_metadata,
                                 control = list(verbose = FALSE))

adjusted_otu_table <- fit_adjust_batch$feature_abd_adj

adj_physeq <- metaphlanToPhyloseq(adjusted_otu_table,
                                  metadat = merged_metadata,
                                  simplenames = FALSE)

random_tree <- rtree(ntaxa(adj_physeq), rooted = T, 
                     tip.label = taxa_names(adj_physeq))

adj_physeq <- merge_phyloseq(adj_physeq, random_tree)


## ----afterMMUPhin-----------------------------------------------------------------------------------------
adj_physeq.clr <- microbiome::transform(adj_physeq, "clr")
adj_aitchison_dist <- distance(adj_physeq.clr, method = "euclidean")
adj_aitchison_pcoa  <- ordinate(adj_physeq.clr, "PCoA", distance = adj_aitchison_dist)

clr.plot.after <- plot_ordination(adj_physeq.clr,
                                  adj_aitchison_pcoa, 
                                  color = "study_name") + 
  geom_point() +
  theme_bw()

clr.plot.after + stat_ellipse(aes(group = study_name), linetype = 2) 



p1 <- physeq.fil %>%
  tax_transform("clr", rank = "Species") %>%
  dist_calc("euclidean") %>%
  ord_calc(method = "PCoA") %>% 
  ord_plot(color = "study_name", 
           shape = "study_name") +
  scale_colour_brewer(palette = "Set2", 
                      aesthetics = c("fill", "colour"), name = "study_name") +
  theme_bw() +
  ggside::geom_xsidedensity(aes(fill = study_name), 
                            alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = study_name), 
                            alpha = 0.5, study_name = FALSE) +
  ggside::theme_ggside_void() +
  ggtitle("Before Correction")

p2 <- adj_physeq %>%
  tax_transform("clr", rank = "Species") %>%
  dist_calc("euclidean") %>%
  ord_calc(method = "PCoA") %>% 
  ord_plot(color = "study_name", 
           shape = "study_name") +
  scale_colour_brewer(palette = "Set2", 
                      aesthetics = c("fill", "colour"), name = "study_name") +
  theme_bw() +
  ggside::geom_xsidedensity(aes(fill = study_name), 
                            alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = study_name), 
                            alpha = 0.5, show.legend = FALSE) +
  ggside::theme_ggside_void() +
  ggtitle("After Correction")

ggarrange(p1, p2, nrow = 1, common.legend = T)


## ----aitchison_after_adonis-------------------------------------------------------------------------------
adj_aitchison_adonis <- vegan::adonis(adj_aitchison_dist ~ study_name,
                                        data = data.frame(sample_data(adj_physeq.clr)))
adj_aitchison_adonis$aov.tab 


## ----save_adj_physeq--------------------------------------------------------------------------------------
saveRDS(adj_physeq, "../data/physeqMHO.rds")
