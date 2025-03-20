
rm(list = ls())

library(curatedMetagenomicData)
library(tidyverse)
library(janitor)
library(readxl)
library(ComplexHeatmap)
library(gtsummary)
library(RColorBrewer)

library(phyloseq)
library(ape)
library(mia)
library(microbiome)

library(MMUPHin)

metadata <- read_csv("metadata_curatedAI4Food_procesado_20240220.csv") %>%
  mutate(subject_id = if_else(study_name == "JieZ_2017", sample_id, subject_id)) %>%
  filter(region == "Europe") %>%
  select(-c(region, western))

df <- sampleMetadata[ ,!colnames(sampleMetadata) %in% colnames(metadata)] %>%
  cbind(sampleMetadata[, c("sample_id", "subject_id")]) %>%
  dplyr::inner_join(metadata, by = c("subject_id", "sample_id"))


tse_ctd <- df  %>%
  returnSamples("relative_abundance", counts = FALSE)
physeq_ctd <- makePhyloseqFromTreeSE(tse_ctd,
                                     assay.type = "relative_abundance")
rm(tse_ctd)


physeq_a4f <- readRDS("metaphlan/physeq_mpa3_MHO_mpa30_20240829.rds")


metaphlanToPhyloseq <- source('./metaphlanToPhyloseq.R')$value
merged_otu_table <- t(dada2::mergeSequenceTables(t(physeq_a4f@otu_table), 
                                                 t(physeq_ctd@otu_table)))


merged_phyloseq <- metaphlanToPhyloseq(merged_otu_table,
                                       simplenames = FALSE)

random_tree <- rtree(ntaxa(merged_phyloseq), rooted = T, 
                     tip.label = taxa_names(merged_phyloseq))

merged_phyloseq <- merge_phyloseq(merged_phyloseq, random_tree)
physeq <- tax_glom(merged_phyloseq, "Species")


sample_data(physeq) <- metadata %>% 
  mutate(isCurated = if_else(study_name == "AI4Food", "No", "Yes")) %>%
  column_to_rownames(var = "sample_id")

physeq

saveRDS(physeq, "mhoCuratedEUR_AI4Food_mpa30_20240829.rds")
