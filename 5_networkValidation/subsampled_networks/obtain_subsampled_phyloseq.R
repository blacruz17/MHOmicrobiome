library(SpiecEasi)
library(phyloseq)
library(genefilter)
library(pulsar)
library(igraph)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)

rm(list = ls())

# ---- Get phyloseq objects for each phenotype ---------------------------------
physeq <- readRDS("physeq_adj.rds")

(physeq.mho  <- subset_samples(physeq, MetObesity == "MHO")) # 19 samples 
(physeq.mhno <- subset_samples(physeq, MetObesity == "MHNO")) # 103 samples 
(physeq.muo  <- subset_samples(physeq, MetObesity == "MUO")) # 432 samples 
(physeq.muno <- subset_samples(physeq, MetObesity == "MUNO")) # 377 samples 


# Random sampling ##############################################################
set.seed(505)

n_iters <- 100
sample_size <- 100

samplingFun <- function(physeqObj, name, sampleSize, n_times){
  
  samples <- lapply(seq_len(n_times),
                    function(x) {
                      data.frame(sample_data(physeqObj)) %>%
                        slice_sample(n = sampleSize, replace = FALSE) %>%
                        rownames()
                    })
  
  
  names(samples) <- paste0(name, "_", 1:sampleSize)
  samplingDF <- data.frame(samples)
  
  return(samplingDF)
}

mhno.samples <- samplingFun(physeq.mhno, "MHNO", sample_size, n_iters)
muno.samples <- samplingFun(physeq.muno, "MUNO", sample_size, n_iters)
muo.samples <- samplingFun(physeq.muo, "MUO", sample_size, n_iters)

sampling.df <- cbind(mhno.samples, muno.samples, muo.samples)
write_csv(sampling.df, "./randomSubsamples.csv")

# Obtain + export phyloseq objects for each random subsample ###################

# MHNO
physeq_list.mhno <- lapply(colnames(mhno.samples), function(cn){
  col <- mhno.samples[, cn]
  pseq.fil <- prune_samples(col, physeq.mhno)
  pseq.fil <- filter_taxa(pseq.fil,
                          flist = filterfun(kOverA(k = 10, A = 0.001)),
                          prune = TRUE)
  
  write_rds(pseq.fil, paste0("./subsampling_100iters/", cn, ".rds"))
})
names(physeq_list.mhno) <- names(mhno.samples)

# MUNO
physeq_list.muno <- lapply(colnames(muno.samples), function(cn){
  col <- muno.samples[, cn]
  pseq.fil <- prune_samples(col, physeq.muno)
  pseq.fil <- filter_taxa(pseq.fil,
                          flist = filterfun(kOverA(k = 10, A = 0.001)),
                          prune = TRUE)
  
  write_rds(pseq.fil, paste0("./subsampling_100iters/", cn, ".rds"))
})
names(physeq_list.muno) <- names(muno.samples)

# MUO
physeq_list.muo <- lapply(colnames(muo.samples), function(cn){
  col <- muo.samples[, cn]
  pseq.fil <- prune_samples(col, physeq.muo)
  pseq.fil <- filter_taxa(pseq.fil,
                          flist = filterfun(kOverA(k = 10, A = 0.001)),
                          prune = TRUE)
  
  write_rds(pseq.fil, paste0("./subsampling_100iters/", cn, ".rds"))
})
names(physeq_list.muo) <- names(muo.samples)

# Export phyloseq taxonomic abundances #########################################
exportMeta <- function(physeq, filename){
  
  # get taxonomic info:
  taxa <- data.frame(physeq@tax_table) %>%
    rownames_to_column(var = "Label")
  
  # get mean abundances:
  mean_ab <- data.frame(
    "Abundance" = apply(data.frame(physeq@otu_table), 1, mean)) %>%
    rownames_to_column(var = "Label")
  
  # join all of them together:
  metadata <- taxa %>%
    left_join(mean_ab, by = "Label") 
  
  # WRITE CSV FILE:
  write_csv(metadata,
            paste0(filename, "_metadata.csv"),
            quote = "needed")
}

lapply(names(physeq_list.mhno), function(pseq_name){
  filename_pseq <- paste0("./subsampling_100iters/", pseq_name)
  exportMeta(physeq_list.mhno[[pseq_name]], filename_pseq)
})

lapply(names(physeq_list.muno), function(pseq_name){
  filename_pseq <- paste0("./subsampling_100iters/", pseq_name)
  exportMeta(physeq_list.muno[[pseq_name]], filename_pseq)
})

lapply(names(physeq_list.muo), function(pseq_name){
  filename_pseq <- paste0("./subsampling_100iters/", pseq_name)
  exportMeta(physeq_list.muo[[pseq_name]], filename_pseq)
  print(filename_pseq)
})
