library(tidyverse)
library(phyloseq)
library(microbiome)
library(genefilter)
library(NetCoMi)

# Data Import & Preprocessing #########
(physeq <- readRDS("../data/physeqMHO.rds")) # 355 taxa
any(rowSums(physeq@otu_table) == 0) 

(physeq.mho <- subset_samples(physeq, MetObesity == "MHO")) # 22 samples
(physeq.muo <- subset_samples(physeq, MetObesity == "MUO")) # 437 samples
(physeq.mhno <- subset_samples(physeq, MetObesity == "MHNO")) # 115 samples
(physeq.muno <- subset_samples(physeq, MetObesity == "MUNO")) # 385 samples

# remove 0s:
(physeq.mho <- filter_taxa(physeq.mho,
                           flist =  filterfun(kOverA(k = 1, A = 0)),
                           prune = TRUE)) # 147 taxa x 14 samples
(physeq.muo <- filter_taxa(physeq.muo,
                           flist =  filterfun(kOverA(k = 1, A = 0)),
                           prune = TRUE)) # 149 taxa x 31 samples
(physeq.mhno <- filter_taxa(physeq.mhno,
                            flist =  filterfun(kOverA(k = 1, A = 0)),
                            prune = TRUE)) # 165 taxa x 24 samples
(physeq.muno <- filter_taxa(physeq.muno,
                            flist =  filterfun(kOverA(k = 1, A = 0)),
                            prune = TRUE)) # 148 taxa x 29 samples

# check:
any(rowSums(physeq.mho@otu_table) == 0)
any(rowSums(physeq.muo@otu_table) == 0)
any(rowSums(physeq.mhno@otu_table) == 0)
any(rowSums(physeq.muno@otu_table) == 0)

library(doParallel)
num_cores <- 4
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# MUO
net.muo <- netConstruct(data = physeq.muo,
                        dataType = "counts",
                        measure = "spieceasi",
                        sparsMethod = "none",
                        nboot = 1000,
                        cores = 4L)
props.muo <- netAnalyze(net.muo,
                        centrLCC = TRUE,
                        clustMethod = "cluster_fast_greedy",
                        hubPar = c("degree", "betweenness"),
                        hubQuant = 0.9,
                        weightDeg = FALSE, normDeg = FALSE)

# MHNO
net.mhno <- netConstruct(data = physeq.mhno,
                         dataType = "counts",
                         measure = "spieceasi",
                         sparsMethod = "none",
                         nboot = 1000,
                         cores = 4L)
props.mhno <- netAnalyze(net.mhno,
                         centrLCC = TRUE,
                         clustMethod = "cluster_fast_greedy",
                         hubPar = c("degree", "betweenness"),
                         hubQuant = 0.9,
                         weightDeg = FALSE, normDeg = FALSE)

# MUNO
net.muno <- netConstruct(data = physeq.muno,
                         dataType = "counts",
                         measure = "spieceasi",
                         sparsMethod = "none",
                         nboot = 1000,
                         cores = 4L)
props.muno <- netAnalyze(net.muno,
                         centrLCC = TRUE,
                         clustMethod = "cluster_fast_greedy",
                         hubPar = c("degree", "betweenness"),
                         hubQuant = 0.9,
                         weightDeg = FALSE, normDeg = FALSE)

# MHO
net.mho <- netConstruct(data = physeq.mho,
                        dataType = "counts",
                        measure = "spieceasi",
                        sparsMethod = "none",
                        nboot = 1000,
                        cores = 4L)
props.mho <- netAnalyze(net.mho,
                        centrLCC = TRUE,
                        clustMethod = "cluster_fast_greedy",
                        hubPar = c("degree", "betweenness"),
                        hubQuant = 0.9,
                        weightDeg = FALSE, normDeg = FALSE)

# Export #######################################################################
# - Edgelists
# - Node Metadata including hubs!

exportNet <- function(net, physeq, props, filename){
  # EDGES
  # Create edge object from the edge list exported by netConstruct()
  edges <- net$edgelist1 %>%
    select(v1, v2, asso) %>%
    rename(Source = v1,
           Target = v2,
           Weight = asso) %>%
    mutate(Type = "Undirected",
           Sign = if_else(Weight < 0, "Negative", "Positive"))
  
  # NODES
  # get taxonomic info:
  taxa <- data.frame(physeq@tax_table) %>%
    rownames_to_column(var = "Label")
  
  # get mean abundances:
  mean_ab <- data.frame(
    "Abundance" = apply(data.frame(physeq@otu_table), 1, mean)) %>%
    rownames_to_column(var = "Label")
  
  # get clusters:
  clusters <- data.frame("Cluster" =  props$clustering$clust1) %>%
    rownames_to_column(var = "Label")
  
  # get hubs:
  hubs <- data.frame(isHub = rep(1, length(props$hubs$hubs1)),
                     "Label" = props$hubs$hubs1)
  
  # join all of them together:
  metadata <- taxa %>%
    left_join(mean_ab, by = "Label") %>%
    left_join(clusters, by = "Label") %>%
    left_join(hubs, by = "Label") %>%
    mutate(isHub = if_else(is.na(isHub), 0, 1))
  
  
  # WRITE CSV FILES:
  write_csv(edges,
            file = paste0(filename, "_edges.csv"),
            quote = "needed")
  write_csv(metadata,
            paste0(filename, "_metadata.csv"),
            quote = "needed")
}

exportNet(net.mho, physeq.mho, props.mho, "../results/mho")
exportNet(net.mhno, physeq.mhno, props.mhno, "../results/mhno")
exportNet(net.muo, physeq.muo, props.muo, "../results/muo")
exportNet(net.muno, physeq.muno, props.muno, "../results/muno")

stopCluster(cl)

# Visualize ####################################################################

#Colors and aesthetics
taxtab.mho <- physeq.mho@tax_table@.Data
phyla.mho <- as.factor(gsub("p__", "", taxtab.mho[, "Phylum"]))

taxtab.mhno <- physeq.mhno@tax_table@.Data
phyla.mhno <- as.factor(gsub("p__", "", taxtab.mhno[, "Phylum"]))

taxtab.muo <- physeq.muo@tax_table@.Data
phyla.muo <- as.factor(gsub("p__", "", taxtab.muo[, "Phylum"]))

taxtab.muno <- physeq.muno@tax_table@.Data
phyla.muno <- as.factor(gsub("p__", "", taxtab.muno[, "Phylum"]))

pal <- c("#264653", "#2A9D8F", "#703d57", "#edafb8", "#e6af2e")
newpal <- colorRampPalette(pal)
pal_12 <- newpal(13)[c(1:5, 7:13)]

levels(phyla.mhno) == levels(phyla.muno)
levels(phyla.mhno) == levels(phyla.muo)
levels(phyla.mhno) == levels(phyla.mho)

pal_mho <- pal_12[c(1:6, 8, 9, 11, 12)]


png("../figures/mho.png",
    width = 16,
    height = 16,
    units = "cm",
    res = 1200)
plot(mho.props, 
     nodeColor = "feature", 
     featVecCol = phyla.mho, 
     colorVec =  pal_mho,
     posCol = pal[2],
     negCol = pal[3],
     borderCol = "gray20",
     nodeSize = "mclr", 
     nodeSizeSpread = 3,
     edgeTranspLow = 30, 
     edgeTranspHigh = 10,
     mar = rep(2, 4),
     repulsion = 1,
     rmSingles = "none",
     highlightHubs = T,
     hubTransp = .7, 
     hubBorderWidth = 3,
     cexLabels=0)
dev.off()


png("../figures/mhno.png",
    width = 16,
    height = 16,
    units = "cm",
    res = 1200)
plot(mhno.props, 
     nodeColor = "feature", 
     featVecCol = phyla.mhno, 
     colorVec =  pal_12,
     posCol = pal[2],
     negCol = pal[3],
     borderCol = "gray20",
     nodeSize = "mclr", 
     nodeSizeSpread = 3,
     edgeTranspLow = 30, 
     edgeTranspHigh = 10,
     mar = rep(2, 4),
     repulsion = 1,
     rmSingles = "none",
     highlightHubs = T,
     hubTransp = .7, 
     hubBorderWidth = 3,
     cexLabels=0)
dev.off()

png("../figures/muno.png",,
    width = 16,
    height = 16,
    units = "cm",
    res = 1200)
plot(muno.props, 
     nodeColor = "feature", 
     featVecCol = phyla.muno, 
     colorVec = pal_12,
     posCol = pal[2],
     negCol = pal[3],
     borderCol = "gray20",
     nodeSize = "mclr", 
     nodeSizeSpread = 3,
     edgeTranspLow = 30, 
     edgeTranspHigh = 10,
     mar = rep(2, 4),
     repulsion = 1,
     rmSingles = "none",
     highlightHubs = T,
     hubTransp = .7, 
     hubBorderWidth = 3,
     cexLabels=0)
dev.off()

png("../figures/muo.png",
    width = 16,
    height = 16,
    units = "cm",
    res = 1200)

plot(muo.props, 
     nodeColor = "feature", 
     featVecCol = phyla.muo, 
     colorVec = pal_12,
     posCol = pal[2],
     negCol = pal[3],
     borderCol = "gray20",
     nodeSize = "mclr", 
     nodeSizeSpread = 3,
     edgeTranspLow = 30, 
     edgeTranspHigh = 10,
     mar = rep(2, 4),
     repulsion = 1,
     rmSingles = "none",
     highlightHubs = T,
     hubTransp = .7, 
     hubBorderWidth = 3,
     cexLabels=0)
dev.off()


