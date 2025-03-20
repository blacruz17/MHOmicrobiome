library(phyloseq)
library(microbiome)
library(ape)


# metadata -------
setwd(dataDir)
metadata_v1 <- read_tsv("metadata_v1.tsv")
metadata_v3 <- read_tsv("metadata_v3.tsv")
metadata_mho <- read_csv(paste0(resultDir, "/tabla_mho_ai4food.csv"))
metadata <- bind_rows(metadata_v1, metadata_v3)
metaphlanToPhyloseq <- source('metaphlanToPhyloseq.R')$value

# create phyloseq object --------------------------
setwd(resultDir)

data <- read_tsv('mpa_v30/profiled_MPA3/merged_abundance_table.txt', skip = 1)
data <- column_to_rownames(data, var = "clade_name")

ids_vol_visita <- data.frame(
  "id_voluntario" = gsub("profiled_V[13]_", "", colnames(data)),
  "visita" = gsub("profiled_V", "", colnames(data)))
ids_vol_visita$visita <- gsub("_[0-9]*", "", ids_vol_visita$visita)

metadata <- inner_join(ids_vol_visita %>%
                         mutate(id_voluntario = as.character(id_voluntario),
                                visita = as.character(visita)),
                       metadata %>%
                         mutate(id_voluntario = as.character(id_voluntario),
                                visita = as.character(visita)),
                       by = c("id_voluntario", "visita")) %>%
  left_join(metadata_mho %>%
              mutate(id_voluntario = as.character(id_voluntario),
                     visita = as.character(visita)),
            by = c("id_voluntario", "visita"))

colnames(data) <- as.numeric(gsub("profiled_V[13]_", "", colnames(data)))
metadata <- column_to_rownames(metadata, var = "id_voluntario")

phyloseqin <- metaphlanToPhyloseq(data, metadata, simplenames = FALSE)

random_tree <- rtree(ntaxa(phyloseqin), rooted = T, 
                     tip.label = taxa_names(phyloseqin))

physeq <- merge_phyloseq(phyloseqin, random_tree)
# tax_glom:
physeq <- tax_glom(physeq, "Species")
# sanity check
summary(colSums(as.matrix(otu_table(physeq)))) # OK!

saveRDS(physeq, 'metaphlan/physeq_mpa3_MHO_mpa30_20240829.rds')
