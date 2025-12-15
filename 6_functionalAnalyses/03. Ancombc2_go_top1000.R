!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Differential abundance analysis of HUMAnN gene families mapped to GO terms
#
# This script analyses HUMAnN gene-family outputs (GO-mapped; relative abundances)
# using ANCOM-BC2. The GO terms are restricted to a candidate set obtained from a
# prior GO-DAG filtering step (02. GO_term_refinement).
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(tibble)
  library(ANCOMBC)
})

# 1) Data import and column name harmonisation
input_data <- read_tsv("Data/genefamilies_merged_unstrat_to_go.tsv")

# Remove technical HUMAnN suffixes from sample column names 
colnames(input_data) <- colnames(input_data) %>%
  str_replace("_hostclean.*Abundance-RPKs$", "") %>%
  str_trim()

# 2) Retain only the samples required for the analysis

id <- readLines("Data/muestraanalisis.txt")

# Retain the gene-family identifier column plus the selected sample IDs.
cols_to_keep <- c("# Gene Family", id)
input_data_filtered <- input_data[, colnames(input_data) %in% cols_to_keep]

# Standardise sample identifiers 
colnames(input_data_filtered) <- colnames(input_data_filtered) %>%
  str_replace("^V[0-9]+_", "") %>%
  str_replace("(_.*)$", "") %>%
  str_trim()


# 3) Restrict rows to GO candidates (derived in the previous notebook/script)

go_candidates <- read_excel(
  "C:/Users/alba.perez/OneDrive - FUNDACION IMDEA-ALIMENTACION/paper_blanca/Data/my_GO_leaves_with_immediate_and_global_ancestor.xlsx"
)

go_candidates <- go_candidates$GO_ID
length(go_candidates)

# Retain candidate GO terms plus HUMAnN bookkeeping rows for context.
input_data_filtered <- input_data_filtered %>%
  filter(`# Gene Family` %in% go_candidates | `# Gene Family` %in% c("UNMAPPED", "UNGROUPED"))

# 4) Summarise UNMAPPED / UNGROUPED proportions across samples

special_rows <- input_data_filtered %>%
  filter(`# Gene Family` %in% c("UNMAPPED", "UNGROUPED"))

special_long <- special_rows %>%
  pivot_longer(cols = -`# Gene Family`, names_to = "sample", values_to = "abundance")

summary_global <- special_long %>%
  group_by(`# Gene Family`) %>%
  summarise(mean_abundance = mean(abundance) * 100)

print(summary_global)

# 5) Exploratory distribution plots for mean abundances across GO terms

pathways_long <- input_data_filtered %>%
  filter(!(`# Gene Family` %in% c("UNMAPPED", "UNGROUPED"))) %>%
  pivot_longer(
    cols = -`# Gene Family`,
    names_to = "sample",
    values_to = "abundance"
  )

rank_curve <- pathways_long %>%
  group_by(`# Gene Family`) %>%
  summarise(mean_abundance = mean(abundance)) %>%
  arrange(desc(mean_abundance)) %>%
  mutate(rank = row_number())

p_rank <- ggplot(rank_curve, aes(x = rank, y = mean_abundance)) +
  geom_line(color = "darkred") +
  scale_y_log10() +
  labs(
    x = "GO term rank (by mean abundance)",
    y = "Mean abundance (log10 scale)"
  )

print(p_rank)

ggsave(
  filename = "Figures/rank_curve.png",
  plot     = p_rank,
  dpi      = 300,
  width    = 12,
  height   = 8,
  units    = "cm",
  bg       = "white"
)

p_hist <- ggplot(rank_curve, aes(x = mean_abundance)) +
  geom_histogram(bins = 50, fill = "darkred", alpha = 0.7) +
  scale_x_log10() +
  labs(x = "Mean abundance (log10 scale)", y = "Frequency")

print(p_hist)

ggsave(
  filename = "Figures/mean_abundance_hist.png",
  plot     = p_hist,
  dpi      = 300,
  width    = 12,
  height   = 8,
  units    = "cm",
  bg       = "white"
)

# Inspect the abundance at rank 1,000 (used as the cut-off below).
rank_1000 <- rank_curve %>% filter(rank == 1000)
print(rank_1000)

# 6) Retain the top 1,000 most abundant GO terms

# Variable names are preserved from the original notebook.
top1000_ids <- rank_curve %>%
  slice_head(n = 1000) %>%
  pull(`# Gene Family`)

input_top1000 <- input_data_filtered %>%
  filter(`# Gene Family` %in% top1000_ids)

# Remove bookkeeping rows prior to modelling.
input_top1000 <- input_top1000 %>%
  filter(!(`# Gene Family` %in% c("UNMAPPED", "UNGROUPED")))

# Set GO IDs as row names; columns remain samples.
input_top1000 <- input_top1000 %>%
  column_to_rownames(var = "# Gene Family")


# 7) Load and align metadata

meta_data <- read_csv("Data/metadatos_20250915.csv", show_col_types = FALSE)

samples_top1000 <- colnames(input_top1000)

meta_data_top1000 <- meta_data %>%
  filter(id_sample %in% samples_top1000) %>%
  arrange(factor(id_sample, levels = samples_top1000)) %>%
  column_to_rownames(var = "id_sample")

# Confirm that sample ordering is identical between the abundance matrix and metadata.
print(identical(colnames(input_top1000), rownames(meta_data_top1000)))

input_data <- input_top1000
input_metadata <- meta_data_top1000

print(all(colnames(input_data) == rownames(input_metadata)))

# 8) ANCOM-BC2 differential abundance analysis

# ANCOM-BC2 expects counts. Here, relative abundances are converted to pseudo-counts.
# This is a pragmatic transformation; interpretation should reflect this choice.
input_counts <- round(1e6 * input_data)

# Set reference levels explicitly to ensure consistent contrasts.
input_metadata$MetObesity <- relevel(factor(input_metadata$MetObesity), ref = "MHNO")
input_metadata$study_name <- relevel(factor(input_metadata$study_name), ref = "ai4food")
input_metadata$sex <- factor(input_metadata$sex)

set.seed(123)

out <- ancombc2(
  data          = input_counts,
  taxa_are_rows = TRUE,  # rows = features (GO terms), columns = samples
  meta_data     = input_metadata,
  fix_formula   = "MetObesity + study_name+ sex + age + bmi_kg_m2",
  group         = "MetObesity",
  p_adj_method  = "holm",
  prv_cut       = 0,
  lib_cut       = 1000,
  s0_perc       = 0.05,
  struc_zero    = TRUE,
  neg_lb        = TRUE,
  alpha         = 0.05,
  global        = FALSE,   # global test for MetObesity
  pairwise      = TRUE,   # all pairwise comparisons
  dunnet        = FALSE,   # each group vs the reference (MHNO)
  n_cl          = 1,
  verbose       = TRUE
)

out <- readRDS("Ancom-bc/diffAbundance_sep25_GO/ANCOMBC2_20251013_allvariable_top1000.RDS")






