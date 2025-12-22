## ----setup, include=FALSE----------------------------------------
knitr::opts_chunk$set(warning = FALSE, echo = TRUE)
rm(list = ls())

# CRAN
library(tidyverse)
library(vegan)
library(gridExtra)
library(ggpubr)
library(patchwork)

# Bioconductor
library(phyloseq)
library(microbiome)
library(genefilter)
library(microViz)
library(MMUPHin)
library(ape)




## ----getPhyloseq-------------------------------------------------
physeq.raw <- readRDS("physeq.rds")
physeq.5m  <- readRDS("physeq_5m.rds")
physeq.8m  <- readRDS("physeq_8m.rds")

metaphlanToPhyloseq <- source('metaphlanToPhyloseq.R')$value



## ----getAlpha----------------------------------------------------
plotAlpha <- function(physeqObj){
  
  # Get Alpha Diversity + Study name
  alpha.res <- alpha(physeqObj, 
                     index = c("chao1")) %>%
                  rownames_to_column(var = "id")
  
  alpha.df <- sample_data(physeqObj) %>%
                data.frame() %>%
                rownames_to_column(var = "id") %>%
                left_join(alpha.res, by = "id")
  
  # Make plots
  alphaPlot <- alpha.df %>%
    ggplot(aes(x = study_name, y = chao1)) +
    geom_boxplot(aes(fill = study_name), 
                 alpha=.6,
                 outlier.shape = NA,
                 outlier.size = 2,
                 color = "black") +
    geom_jitter(aes(shape = study_name, color = study_name),
                width = .2,
                size = .8) +
    theme_bw()
  
  return(alphaPlot)
}

alpha.p.raw <- plotAlpha(physeq.raw) + ggtitle("All reads")
alpha.p.5m  <- plotAlpha(physeq.5m)  + ggtitle("5M reads")
alpha.p.8m  <- plotAlpha(physeq.8m)  + ggtitle("8M reads")

ggarrange(alpha.p.raw, alpha.p.5m, alpha.p.8m,
          nrow = 1,
          common.legend = TRUE,
          legend = "bottom")


## ----functionFilters---------------------------------------------
select_filter_phyloseq <- function(physeq,
                                   cutoff_abundance,
                                   cutoff_prevalence) {
  
  # Param combinations
  tb_param <- tidyr::crossing(cutoff_abundance, cutoff_prevalence)
  
  # 1. Filtering using all possible combinations
  tb_filter <- purrr::map_dfr(seq_len(nrow(tb_param)), function(i_param) {
    cutoff_abundance <- tb_param$cutoff_abundance[i_param]
    cutoff_prevalence <- tb_param$cutoff_prevalence[i_param]
    
    tb_presence <- physeq %>%
      group_by(Species, study_name) %>%
      summarise(
        n = n(),
        present = sum(Abundance > cutoff_abundance) > n * cutoff_prevalence,
        .groups = "drop"
      )
    
    tb_presence %>%
      group_by(Species) %>%
      summarise(n_study = sum(present), .groups = "drop") %>%
      mutate(
        cutoff_abundance = cutoff_abundance,
        cutoff_prevalence = cutoff_prevalence,
        feature_level = "Species"
      )
  })
  
  # 2. Summary table
  tb_summary <- tb_filter %>%
    group_by(feature_level, cutoff_abundance, cutoff_prevalence) %>%
    summarise(
      nfeature_nonUnique = sum(n_study > 1),
      nfeature_unique = sum(n_study == 1),
      nfeature_total = nfeature_nonUnique + nfeature_unique,
      .groups = "drop"
    ) %>%
    arrange(feature_level, desc(cutoff_abundance), desc(cutoff_prevalence))
  
  # 3. Select best filter
  best_filter <- tb_summary %>%
    arrange(desc(nfeature_nonUnique)) %>%
    slice(1)
  
  return(list(
    tb_filter = tb_filter,
    tb_summary = tb_summary,
    best_filter = best_filter
  ))
}


## ----psmelt------------------------------------------------------
physeq.raw.mt <- psmelt(physeq.raw)
physeq.5m.mt  <- psmelt(physeq.5m)
physeq.8m.mt  <- psmelt(physeq.8m)


## ----tryfilters--------------------------------------------------
cutoff_abundance <- c(1e-5, 5e-5, 
                      1e-4, 5e-4,
                      1e-3, 5e-3, 
                      1e-2, 5e-2)
cutoff_prevalence <- c(0.05, 0.1)

fil.raw <- select_filter_phyloseq(physeq.raw.mt, 
                                  cutoff_abundance = cutoff_abundance,
                                  cutoff_prevalence = cutoff_prevalence)
fil.raw$tb_summary; fil.raw$best_filter

fil.5m  <- select_filter_phyloseq(physeq.5m.mt, 
                                  cutoff_abundance = cutoff_abundance,
                                  cutoff_prevalence = cutoff_prevalence)
fil.5m$tb_summary; fil.5m$best_filter

fil.8m  <- select_filter_phyloseq(physeq.8m.mt, 
                                  cutoff_abundance = cutoff_abundance,
                                  cutoff_prevalence = cutoff_prevalence)
fil.8m$tb_summary; fil.8m$best_filter


## ----Filter------------------------------------------------------
filter_physeq <- function(physeq, physeq.fil.info){
  
  sp.keep <- physeq.fil.info$tb_filter %>% 
    dplyr::filter(cutoff_abundance == physeq.fil.info$best_filter$cutoff_abundance,
                  cutoff_prevalence == physeq.fil.info$best_filter$cutoff_prevalence,
                  n_study >= 1) %>% 
    dplyr::pull(Species) 
  
  otu.keep <- data.frame(tax_table(physeq)) %>%
    filter(Species %in% sp.keep) %>%
    rownames()
  
  
  physeq.filtered <-  phyloseq::prune_taxa(otu.keep, physeq) 
  return(physeq.filtered)
}

physeq.raw.fil <- filter_physeq(physeq.raw, fil.raw)
physeq.5m.fil  <- filter_physeq(physeq.5m, fil.5m)
physeq.8m.fil  <- filter_physeq(physeq.8m, fil.8m)


## ----clrTrans----------------------------------------------------
physeq.raw.clr <- microbiome::transform(physeq.raw.fil, "clr")
physeq.5m.clr  <- microbiome::transform(physeq.5m.fil, "clr")
physeq.8m.clr  <- microbiome::transform(physeq.8m.fil, "clr")


## ----atichi_study------------------------------------------------
# No rarefaction
aitchison_dist  <- distance(physeq.raw.clr, method = "euclidean")
aitchison_pcoa  <- ordinate(physeq.raw.clr, "PCoA", distance = aitchison_dist)
# 5 million
aitchison_dist_5m  <- distance(physeq.5m.clr, method = "euclidean")
aitchison_pcoa_5m  <- ordinate(physeq.5m.clr, "PCoA", distance = aitchison_dist_5m)
# 8 million
aitchison_dist_8m  <- distance(physeq.8m.clr, method = "euclidean")
aitchison_pcoa_8m  <- ordinate(physeq.8m.clr, "PCoA", distance = aitchison_dist_8m)


## ----pcoa_study--------------------------------------------------
plot_ordination(physeq.raw.clr, aitchison_pcoa, color = "study_name") +
  geom_point() + theme_bw() + ggtitle("All reads")

plot_ordination(physeq.5m.clr, aitchison_pcoa_5m, color = "study_name") +
  geom_point() + theme_bw() + ggtitle("5 mill reads")

plot_ordination(physeq.8m.clr, aitchison_pcoa_8m, color = "study_name") +
  geom_point() + theme_bw() + ggtitle("8 mill reads")


## ----AdonisStudyName---------------------------------------------
aitchison_adonis_study <- vegan::adonis2(
  aitchison_dist ~ study_name,
  data = data.frame(sample_data(physeq.raw.clr)))
aitchison_adonis_study 


## ----AdonisStudyName_5m------------------------------------------
aitchison_adonis_study_5m <- vegan::adonis2(
  aitchison_dist_5m ~ study_name,
  data = data.frame(sample_data(physeq.5m.clr)))
aitchison_adonis_study_5m 


## ----AdonisStudyName_8m------------------------------------------
aitchison_adonis_study_8m <- vegan::adonis2(
  aitchison_dist_8m ~ study_name,
  data = data.frame(sample_data(physeq.8m.clr)))
aitchison_adonis_study_8m 


## ----mmuphin-----------------------------------------------------
# converts to numbers 0-1:
merged_otu <- otu_table(physeq.8m.fil) %>% 
  apply(2, function(x) x/100)
merged_metadata <- data.frame(sample_data(physeq.8m.fil))

# batch effect adjustment:
fit_adjust_batch <- adjust_batch(feature_abd = merged_otu,
                                 batch = "study_name",
                                 covariates = "MetObesity",
                                 data = merged_metadata,
                                 control = list(verbose = FALSE))

adj_physeq <- metaphlanToPhyloseq(fit_adjust_batch$feature_abd_adj,
                                  metadat = merged_metadata,
                                  simplenames = FALSE)

random_tree <- rtree(ntaxa(adj_physeq), rooted = T, 
                     tip.label = taxa_names(adj_physeq))

adj_physeq <- merge_phyloseq(adj_physeq, random_tree)


## ----plotmmuphin-------------------------------------------------
adj_physeq.clr <- microbiome::transform(adj_physeq, "clr")
adj_aitchison_dist <- distance(adj_physeq.clr, method = "euclidean")
adj_aitchison_pcoa  <- ordinate(adj_physeq.clr, "PCoA", distance = adj_aitchison_dist)

clr.plot.after <- plot_ordination(adj_physeq.clr,
                                  adj_aitchison_pcoa, 
                                  color = "study_name") + 
  geom_point() +
  theme_bw()

clr.plot.after + stat_ellipse(aes(group = study_name), linetype = 2) 


## ----statsmmuphin------------------------------------------------
adj_aitchison_adonis_study <- vegan::adonis2(
  adj_aitchison_dist ~ study_name,
  data = data.frame(sample_data(adj_physeq.clr)))
adj_aitchison_adonis_study 


## ----statsmmuphinMHO---------------------------------------------
adj_aitchison_adonis_mho <- vegan::adonis2(
  adj_aitchison_dist ~ MetObesity,
  data = data.frame(sample_data(adj_physeq.clr)))
adj_aitchison_adonis_mho 


## ----mmuphin5m---------------------------------------------------
# converts to numbers 0-1:
merged_otu <- otu_table(physeq.5m.fil) %>% 
  apply(2, function(x) x/100)
merged_metadata <- data.frame(sample_data(physeq.5m.fil))

# batch effect adjustment:
fit_adjust_batch <- adjust_batch(feature_abd = merged_otu,
                                 batch = "study_name",
                                 covariates = "MetObesity",
                                 data = merged_metadata,
                                 control = list(verbose = FALSE))

adj_physeq <- metaphlanToPhyloseq(fit_adjust_batch$feature_abd_adj,
                                  metadat = merged_metadata,
                                  simplenames = FALSE)

random_tree <- rtree(ntaxa(adj_physeq), rooted = T, 
                     tip.label = taxa_names(adj_physeq))

adj_physeq <- merge_phyloseq(adj_physeq, random_tree)


## ----plotmmuphin5m-----------------------------------------------
adj_physeq.clr <- microbiome::transform(adj_physeq, "clr")
adj_aitchison_dist <- distance(adj_physeq.clr, method = "euclidean")
adj_aitchison_pcoa  <- ordinate(adj_physeq.clr, "PCoA", distance = adj_aitchison_dist)

clr.plot.after <- plot_ordination(adj_physeq.clr,
                                  adj_aitchison_pcoa, 
                                  color = "study_name") + 
  geom_point() +
  theme_bw()

clr.plot.after + stat_ellipse(aes(group = study_name), linetype = 2) 


## ----statsmmuphin5m----------------------------------------------
adj_aitchison_adonis_study <- vegan::adonis2(
  adj_aitchison_dist ~ study_name,
  data = data.frame(sample_data(adj_physeq.clr)))
adj_aitchison_adonis_study 


## ----statsmmuphinMHO5m-------------------------------------------
adj_aitchison_adonis_mho <- vegan::adonis2(
  adj_aitchison_dist ~ MetObesity,
  data = data.frame(sample_data(adj_physeq.clr)))
adj_aitchison_adonis_mho 


## ----mmuphinraw--------------------------------------------------
# converts to numbers 0-1:
merged_otu <- otu_table(physeq.raw.fil) %>% 
  apply(2, function(x) x/100)
merged_metadata <- data.frame(sample_data(physeq.raw.fil))

# batch effect adjustment:
fit_adjust_batch <- adjust_batch(feature_abd = merged_otu,
                                 batch = "study_name",
                                 covariates = "MetObesity",
                                 data = merged_metadata,
                                 control = list(verbose = FALSE))

adj_physeq <- metaphlanToPhyloseq(fit_adjust_batch$feature_abd_adj,
                                  metadat = merged_metadata,
                                  simplenames = FALSE)

random_tree <- rtree(ntaxa(adj_physeq), rooted = T, 
                     tip.label = taxa_names(adj_physeq))

adj_physeq <- merge_phyloseq(adj_physeq, random_tree)


## ----plotmmuphinraw----------------------------------------------
adj_physeq.clr <- microbiome::transform(adj_physeq, "clr")
adj_aitchison_dist <- distance(adj_physeq.clr, method = "euclidean")
adj_aitchison_pcoa  <- ordinate(adj_physeq.clr, "PCoA", distance = adj_aitchison_dist)

clr.plot.after <- plot_ordination(adj_physeq.clr,
                                  adj_aitchison_pcoa, 
                                  color = "study_name") + 
  geom_point() +
  theme_bw()

clr.plot.after + stat_ellipse(aes(group = study_name), linetype = 2) 


## ----statsmmuphinraw---------------------------------------------
adj_aitchison_adonis_study <- vegan::adonis2(
  adj_aitchison_dist ~ study_name,
  data = data.frame(sample_data(adj_physeq.clr)))
adj_aitchison_adonis_study 


## ----statsmmuphinMHOraw------------------------------------------
adj_aitchison_adonis_mho <- vegan::adonis2(
  adj_aitchison_dist ~ MetObesity,
  data = data.frame(sample_data(adj_physeq.clr)))
adj_aitchison_adonis_mho 


## ----savePhyseqRaw-----------------------------------------------
saveRDS(adj_physeq, "physeq_adj.rds")

