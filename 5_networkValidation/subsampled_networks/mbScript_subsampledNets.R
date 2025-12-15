#!/usr/bin/env Rscript

# 0. Load libraries and set wd ################################################

library(SpiecEasi)
library(phyloseq)
library(pulsar)
library(igraph)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)

rm(list = ls())

# 1. Get arguments + import phyloseq ###########################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Use: Rscript process_spiec.R <file.rds>")
}

filename <- args[1]

physeq <- readRDS(filename)

outdir <- paste0(dirname(filename), "./subsampling_100iters")

base_file <- gsub(".rds", "", basename(filename))

# 2. Generate MB SpiecEasi network #############################################
se_mb <- spiec.easi(physeq,
                    method = "mb",
                    sel.criterion = "stars",
                    pulsar.select = TRUE,
                    pulsar.params = list(rep.num = 100,
                                         subsample.ratio = 0.8,
                                         thresh = 0.05,
                                         seed=10010, 
                                         ncores = 10),
                    verbose = TRUE,
                    nlambda = 30)

write_rds(se_mb, 
          paste0(outdir, "/", base_file, "_se_mb.rds"))
          

# 3. Get + export igraph/graphml object ########################################
# get igraph object
ig_mb <- adj2igraph(getRefit(se_mb), 
                    vertex.attr = list(name = taxa_names(physeq)))

# export:
write_graph(ig_mb, 
            paste0(outdir, "/", base_file, "_ig_mb.graphml"),
            format = "graphml")


# 4. Export edge weigths #######################################################
# Edge weights


get_edge_weights <- function(se_obj, ig_obj, physeq_obj){
  
  edge.mx  <- as.matrix(symBeta(getOptBeta(se_obj), mode='maxabs'))
  
  rownames(edge.mx) <- rownames(otu_table(physeq_obj))
  colnames(edge.mx) <- rownames(otu_table(physeq_obj))
  
  se_edgelist <- data.frame(as_edgelist(ig_obj)) %>%
    rename(OTU_1 = X1, OTU_2 = X2)
  
  edge_df <- as.data.frame(as.matrix(edge.mx)) %>% 
    rownames_to_column(var = "OTU_1") %>% 
    reshape2::melt(id.vars = "OTU_1") %>%
    rename(OTU_2 = variable)
  
  edge_weightlist <- se_edgelist %>%
    left_join(edge_df, by = c("OTU_1", "OTU_2"))
  
  return(edge_weightlist)
}

edge_weights <- get_edge_weights(se_mb, ig_mb, physeq)
colnames(edge_weights)[3] <- "weight"

write_csv(edge_weights, 
          paste0(outdir, "/", base_file, "_edge_weights.csv"))


write_csv(edge_support, 
          paste0(outdir, "/", base_file, "_edge_support.csv"))