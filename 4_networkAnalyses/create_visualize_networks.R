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


# ---- Generar objetos phyloseq para cada grupo --------------------------------
physeq <- readRDS("physeq_adj.rds")
physeq

(physeq.mho <- subset_samples(physeq, MetObesity == "MHO")) # 19 samples 
(physeq.mhno <- subset_samples(physeq, MetObesity == "MHNO")) # 103 samples 
(physeq.muo <- subset_samples(physeq, MetObesity == "MUO")) # 432 samples 
(physeq.muno <- subset_samples(physeq, MetObesity == "MUNO")) # 377 samples 


# remove 0s:
(physeq.mho <- filter_taxa(physeq.mho,
                           flist =  filterfun(kOverA(k = 2, A = 0.001)),
                           prune = TRUE)) 
(physeq.mhno <- filter_taxa(physeq.mhno,
                            flist =  filterfun(kOverA(k = 10, A = 0.001)),
                            prune = TRUE)) 
(physeq.muo <- filter_taxa(physeq.muo,
                           flist =  filterfun(kOverA(k = 43, A = 0.001)),
                           prune = TRUE))
(physeq.muno <- filter_taxa(physeq.muno,
                            flist =  filterfun(kOverA(k = 38, A = 0.001)),
                            prune = TRUE))

# export physeq metadata ----------------------
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

exportMeta(physeq.mhno, "./mhno_metadata")
exportMeta(physeq.mho,  "./mho_metadata")
exportMeta(physeq.muno, "./muno_metadata")
exportMeta(physeq.muo,  "./muo_metadata")

# Build spiecEasi networks -----------------------------------------------------

set.seed(42)

mhno_mb <- spiec.easi(physeq.mhno,
                      method = "mb",
                      sel.criterion = "stars",
                      pulsar.select = TRUE,
                      pulsar.params = list(rep.num = 100,
                                           subsample.ratio = 0.8,
                                           thresh = 0.05,
                                           seed=10010, 
                                           ncores = 4
                      ),
                      verbose = TRUE,
                      nlambda = 30)
saveRDS(mhno_mb, "./mb_mhno_net.rds")


ig.mhno.mb <- adj2igraph(getRefit(mhno_mb),  
                         vertex.attr=list(name=taxa_names(physeq.mhno)))
write_graph(ig.mhno.mb, "mhno_mb_net.graphml", format="graphml")


mho_mb <- spiec.easi(physeq.mho,
                     method = "mb",
                     sel.criterion = "stars",
                     pulsar.select = TRUE,
                     pulsar.params = list(rep.num = 100,
                                          subsample.ratio = 0.8,
                                          thresh = 0.05,
                                          seed=10010, 
                                          ncores = 4
                     ),
                     verbose = TRUE,
                     nlambda = 30)
saveRDS(mho_mb, "./mb_mho_net.rds")

ig.mho.mb <- adj2igraph(getRefit(mho_mb),  
                        vertex.attr=list(name=taxa_names(physeq.mho)))
write_graph(ig.mho.mb, "mho_mb_net.graphml", format="graphml")

muno_mb <- spiec.easi(physeq.muno,
                      method = "mb",
                      sel.criterion = "stars",
                      pulsar.select = TRUE,
                      pulsar.params = list(rep.num = 100,
                                           subsample.ratio = 0.8,
                                           thresh = 0.05,
                                           seed=10010, 
                                           ncores = 4
                      ),
                      verbose = TRUE,
                      nlambda = 30)
saveRDS(muno_mb, "./mb_muno_net.rds")

ig.muno.mb <- adj2igraph(getRefit(muno_mb),  
                         vertex.attr=list(name=taxa_names(physeq.muno)))
write_graph(ig.muno.mb, "muno_mb_net.graphml", format="graphml")

muo_mb <- spiec.easi(physeq.muo,
                     method = "mb",
                     sel.criterion = "stars",
                     pulsar.select = TRUE,
                     pulsar.params = list(rep.num = 100,
                                          subsample.ratio = 0.8,
                                          thresh = 0.05,
                                          seed=10010, 
                                          ncores = 4
                     ),
                     verbose = TRUE,
                     nlambda = 30)
saveRDS(muo_mb, "./mb_muo_net.rds")

ig.muo.mb <- adj2igraph(getRefit(muo_mb),  
                        vertex.attr=list(name=taxa_names(physeq.muo)))
write_graph(ig.muo.mb, "muo_mb_net.graphml", format="graphml")

## Edge stability matrices -------------------------------------------------
get_edge_support <- function(se_obj, ig_obj, physeq_obj){
  
  edge_st <- getOptMerge(se_obj)
  rownames(edge_st) <- rownames(otu_table(physeq_obj))
  colnames(edge_st) <- rownames(otu_table(physeq_obj))
  
  se_edgelist <- data.frame(as_edgelist(ig_obj)) %>%
    rename(OTU_1 = X1, OTU_2 = X2)
  edge_stab_df <- as.data.frame(as.matrix(edge_st)) %>% 
    rownames_to_column(var = "OTU_1") %>% 
    reshape2::melt(id.vars = "OTU_1") %>%
    rename(OTU_2 = variable)
  
  edge_support <- se_edgelist %>%
    left_join(edge_stab_df, by = c("OTU_1", "OTU_2"))
  
  return(edge_support)
}

# mb
mhno_support <- get_edge_support(mhno_mb, ig.mhno.mb, physeq.mhno)
mho_support <- get_edge_support(mho_mb, ig.mho.mb, physeq.mho)
muno_support <- get_edge_support(muno_mb, ig.muno.mb, physeq.muno)
muo_support <- get_edge_support(muo_mb, ig.muo.mb, physeq.muo)


supportDF <- bind_rows(mhno_support %>% mutate(Network = "MHNO"),
                       mho_support %>% mutate(Network = "MHO"),
                       muno_support %>% mutate(Network = "MUNO"),
                       muo_support %>% mutate(Network = "MUO"))
pal <- c("#264653", "#2A9D8F", "#703d57", "#edafb8", "#e6af2e")

edge.p <- supportDF %>%
  ggplot(aes(x = value, fill = Network)) +
  geom_histogram(binwidth = .05) +
  scale_fill_manual(values = pal) +
  theme_bw() +
  facet_wrap(~ Network, scales = "free_y") +
  labs(x = "Edge Bootstrap Support (%)",
       y = "Count")
edge.p  
edge.p

ggsave("./edgeSupport.png",
       edge.p,
       bg = "white",
       height = 8,
       width = 12,
       units = "cm",
       dpi = 800)

## stars curves ----------------------------

png("starsCurves.png",
    width = 1600, height = 1600, res = 300) 

par(mfrow = c(2, 2),    
    mar = c(4, 4, 2, 1)) 

plot(mhno_mb$select, legends = FALSE, main = "MHNO")
plot(mho_mb$select,  legends = FALSE, main = "MHO")
plot(muno_mb$select, legends = FALSE, main = "MUNO")
plot(muo_mb$select,  legends = FALSE, main = "MUO")

dev.off()

## Edge weights ------------------------
edgeweights <- function(physeqObj, spiecObj, igraphObj){
  
  weight.mat <- symBeta(getOptBeta(spiecObj), mode = "maxabs")
  colnames(weight.mat) <- taxa_names(physeqObj)
  rownames(weight.mat) <- taxa_names(physeqObj)
  
  edge.list <- ends(igraphObj, E(igraphObj))
  
  edge.weight <- apply(edge.list, 1, function(x){
    weight.mat[x[1], x[2]]
  })
  
  edge.df <- data.frame(
    v1 = edge.list[,1],
    v2 = edge.list[,2],
    edgeweight = edge.weight,
    stringsAsFactors = FALSE
  )
  
  return(edge.df)
}

mhno.ew <- edgeweights(physeq.mhno, mhno_mb, ig.mhno.mb)
mho.ew <- edgeweights(physeq.mho, mho_mb, ig.mho.mb)
muno.ew <- edgeweights(physeq.muno, muno_mb, ig.muno.mb)
muo.ew <- edgeweights(physeq.muo, muo_mb, ig.muo.mb)


write_csv(mhno.ew, "edgeweights_mhno.csv")
write_csv(mho.ew, "edgeweights_mho.csv")
write_csv(muno.ew, "edgeweights_muno.csv")
write_csv(muo.ew, "edgeweights_muo.csv")


# Visualize networks -----------------------------------------------------------

# WARNING: keystone taxa information is necessary to run this section of code.
# run the needed script and then come back!
keystone.df <- read_csv("keystoneTaxa.csv")

# 1 - Obtain consistent color palette.
# we want to get ALL phyla present

all.phyla <- unique(c(
  # all samples
  data.frame(tax_table(mhno.ps))$Phylum,
  data.frame(tax_table(mho.ps))$Phylum,
  data.frame(tax_table(muno.ps))$Phylum,
  data.frame(tax_table(muo.ps))$Phylum,
  
))

# build color palette 

okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00",
  "#CC79A7", "#999999"
)

pal <- colorRampPalette(okabe_ito)(length(all.phyla))


pal.df <- data.frame(cbind(pal, all.phyla))
colnames(pal.df) <- c("Color", "Phylum")
pal.df


## Helper function -----------------
plotNetwork <- function(physeqObj, igraphObj, spiecObj, keystoneDF, netName,
                        seed = 42){
  # color matching
  col.df <- data.frame(V(igraphObj)) %>%
    rownames_to_column(var = "OTU") %>%
    left_join(data.frame(tax_table(physeqObj)) %>%
                rownames_to_column(var = "OTU"),
              by = "OTU") %>%
    left_join(pal.df, by = "Phylum")
  
  # abundances
  ab.df <- if (max(otu_table(physeqObj)) > 1) {
    apply(data.frame(otu_table(physeqObj))/ 100, 1, mean) 
  } else {
    apply(data.frame(otu_table(physeqObj)), 1, mean)
  }
  
  # edge weights
  weight.mat <- symBeta(getOptBeta(spiecObj))
  colnames(weight.mat) <- taxa_names(physeqObj)
  rownames(weight.mat) <- taxa_names(physeqObj)
  
  edge.weights <- sapply(E(igraphObj), function(e) {
    v1 <- ends(igraphObj, e)[1]
    v2 <- ends(igraphObj, e)[2]
    weight.mat[v1, v2]
  })
  
  # keystone taxa
  key.taxa <- keystoneDF %>%
    filter(Network == netName, keystone == TRUE) %>%
    pull(OTU)
  
  (is.keystone <- V(igraphObj)$name %in% key.taxa)
  
  
  vertex.frame.color <- ifelse(is.keystone, "black", "gray20")
  vertex.frame.width <- ifelse(is.keystone, 4, 1)
  
  # adjust transparency
  node.color <- adjustcolor(col.df$Color, alpha.f = 0.6)
  node.color[is.keystone] <- adjustcolor(col.df$Color[is.keystone], alpha.f = 1)
  
  
  # plot network
  set.seed(seed)
  am.coord <- layout_with_fr(igraphObj)
  
  plot(igraphObj, layout=am.coord, 
       vertex.label=NA,
       vertex.color = node.color,
       vertex.size = abs(log(ab.df)),
       vertex.frame.color = vertex.frame.color,
       vertex.frame.width = vertex.frame.width,
       
       edge.color = adjustcolor(ifelse(edge.weights > 0,
                                       "gray20", # + associations
                                       "#E41A1C"  # - associations
       ) , 
       alpha.f = 0.6),
       edge.width = 2
  )
  
}

# Export networks ------------------
exportNetworkPNG <- function(physeqObj, igraphObj, spiecObj, keystoneDF, netName,
                             filename, width = 3000, height = 3000, res = 300, seed = 42) {

  png(filename,
      width = width, height = height, res = res, units = "px")
  par(omi=c(0,0,0,0), mgp=c(0,0,0),mar=c(0,0,0,0))

  plotNetwork(physeqObj, igraphObj, spiecObj, keystoneDF, netName, seed = seed)
  
  dev.off()
}


exportNetworkPNG(physeq.mhno, ig.mhno.mb, mhno_mb, keystone.df, "MHNO", "mhno.png")
exportNetworkPNG(physeq.mho, ig.mho.mb, mho_mb, keystone.df, "MHO", "mho.png")
exportNetworkPNG(physeq.muno, ig.muno.mb, muno_mb, keystone.df, "MUNO", "muno.png")
exportNetworkPNG(physeq.muo, ig.muo.mb, muo_mb, keystone.df, "MUO", "muo.png")