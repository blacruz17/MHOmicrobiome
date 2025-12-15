#!/bin/bash

# ------------------------------------------------------------------------------
# HUMAnN post-processing pipeline (Gene Families → Relative Abundance → GO terms)
#
# This script consolidates per-sample HUMAnN gene family outputs, separates
# stratified/unstratified contributions, normalises the unstratified table to
# relative abundances, and finally regroups UniRef90 gene families to Gene
# Ontology (GO) terms for downstream analysis.
# ------------------------------------------------------------------------------

module load humann/4.0

base_input="/home/proyectos/imdeaalim/alba/paper_blanca/all_results/humann_mho"
base_output="/home/proyectos/imdeaalim/alba/paper_blanca"
split_dir="$base_output/split_tables"
mkdir -p "$split_dir"

# 1) Merge all per-sample genefamilies into one table 
humann_join_tables \
  -i "$base_input" \
  -o "$base_output/genefamilies_merged.tsv" \
  --file_name genefamilies

# 2) Split merged table into unstratified + stratified outputs
humann_split_stratified_table \
  --input "$base_output/genefamilies_merged.tsv" \
  --output "$split_dir"

# 3) Normalize the unstratified table to relative abundances
humann_renorm_table \
  -i "$split_dir/genefamilies_merged_unstratified.tsv" \
  -u relab \
  -o "$base_output/genefamilies_merged_unstrat_relab.tsv"

# 4) Regroup UniRef90 gene families to GO terms
humann_regroup_table \
  -i "$base_output/genefamilies_merged_unstrat_relab.tsv" \
  -g uniref90_go \
  -o "$base_output/genefamilies_merged_unstrat_to_go.tsv"
