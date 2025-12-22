library(phyloseq)
library(readr)

prefix = "muo"
cols <- c("sex", "age", "bmi_kg_m2")
ps <- read_rds("physeq_adj.rds")

sd <- methods::as(phyloseq::sample_data(ps), "data.frame")
na_col <- any(is.na(sd[ , cols, drop = FALSE]))
cat("Na: ", na_col,  "\n")

sd_sel <- as.data.frame(sd[, cols, drop = FALSE])
sd_out <- paste0(prefix, "_meta_data_mpa4.tsv")
write.table(sd_sel, file = sd_out, sep = "\t", quote = FALSE, col.names = NA)

otu <- as(otu_table(ps), "matrix")
otu <- t(otu)
otu_out <- paste0(prefix, "_datampa4.tsv")
write.table(otu, file = otu_out, sep = "\t", quote = FALSE, col.names = NA)


samples_sd  <- rownames(sd_sel)
samples_otu <- rownames(otu)
same_set   <- setequal(samples_sd, samples_otu)
same_order <- identical(samples_sd, samples_otu)
cat("Same sample set: ", same_set,  "\n")
cat("Same column order:   ", same_order, "\n")
