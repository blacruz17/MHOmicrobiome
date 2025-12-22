# MHOmicrobiome
Microbial community analyses on the topic of metabolically healthy obesity.

This repository contains a set of scripts used to preprocess sequencing data, analyze metadata, explore gut microbiome profiles, and perform network-based analyses. Our code focuses on the comparison of four groups encompassing different metabolic health and obesity phenotypes.

## Table of contents
* [1.- FASTQ file preprocessing](#1-fastq-file-preprocessing)
* [2.- Exploratory metadata analyses](#2-exploratory-metadata-analyses)
* [3.- Exploratory gut microbiome analyses](#3-exploratory-gut-microbiome-analyses)
* [4.- Network analyses](#4-network-analyses)
* [5.- Network validation](#5-network-validation)
* [6.- Functional analyses](#6-functional-analyses)


## 1 Fastq file preprocessing

Scripts used to preprocess fastq files. All files were preprocessed using `fastp` for quality control, including filtering of host reads. This was done using a Nextflow script (`preprocess.nf`). Bash scripts for running this pipeline on paired-end or single-end reads are provided.

## 2 Exploratory metadata analyses
R scripts used to generate tables & supplementary figures related to subject metadata:

- `phenotypesA4F_check.R`, `phenotypesMetaCardis_check.R`, `phenotypesKarlsson_Feng_check.R`: scripts used to process original files from publications' supplementary information or from original GitHub / Zenodo repositories and assign MHNO/MHO/MUNO/MUO phenotypes.
- `linearModels_Metadata.R`: compares biochemical and anthropometrical metadata among the four groups, generating a table with all comparisons as well as boxplots.
- `supp_heatmap.R`: explores all datasets within `curatedMetagenomicData` (version 3.8) to look for studies with stool samples from adult subjects not taking antibiotics, and where at least one of our variables of interest to describe obesity and metabolic phenotypes was available. 

Library requirements are specified in the beggining of each R script.

## 3 Exploratory gut microbiome analyses
R scripts for batch effect correction, alpha and beta diversity calculation, and differential abundance analyses:

- `mmuphin_batch_correction.R`: this is the first of the 3 scripts to be executed. It performs batch effect correction caused by the study of origin while controlling for the effect of metabolic health and obesity using MMUPHin (version 1.14). Different prevalence filters are tested.
- `alpha_beta.R`: alpha- and beta-diversity exploration.
- `ancombc2.R`: uses Analysis of Compositions of Microbiomes with Bias Correction 2 (ANCOM-BC2), an extension of the ANCOM-BC methodology that can be implemented in datasets with multiple groups (version 2.2.2). 

Library requirements are specified in the beggining of each R script.

## 4 Network analyses

- `create_visualize_networks.R`: co-occurrence network generation & visualization using SpiecEasi R library. Exports networks and their weights as csv files for their use in subsequent analyses. Edge bootstrap support & StARS model selection curves are explored. _Note:_ network visualization requires information on keystone taxa obtained from `networkMetrics.py`
- `networkMetrics.py`: calculates basic network properties and keystone taxa with NetworkX.
- `attackFramework.py`: performs attacks on networks:
  - random attacks, where 1000 randomizations are performed, using a seed for reproducibility
  - targeted attacks based on degree, betweenness, and mean abundance within the dataset.
- `other_network_plots.R`: generates plots based on output files from previous scripts.

These scripts can be run locally and Python scripts can also be easily executed in cloud services such as Google Colab.

## 5 Network validation
Network validation scripts:
- `flashWeave`: scripts for pre-FlashWeave processing (`preprocess_for_flashweave.R`), FlashWeave network construction (`run_flashweave.jl`) and subsequent analysis (`analyze_fw_graphs.py`)
- `subsampled_networks`: scripts to build subsampled networks and their analysis:
  - `obtain_subsampled_phyloseq.R`: obtain n subsampled phyloseq objects with n samples each. 
  - `mbScript_subsampledNets.R`: reads subsampled phyloseq objects and builds Spiec-Easi networks. This script is run with `run_mbScript_subsampledNets.sh`
  - `metrics_subsampledNets.py`, `attacks_subsampledNets.py`: calculate metrics & execute attacks on sets of subsampled networks
  - `analyze_subsampled_metrics.R`, `keystone_plot_subsamples.R`: analyze results from previous Python scripts

All validation tests performed for AI4Food and MetaCardis cohorts can be run with the scripts provided here by selecting only the samples of interest.

## 6 Functional analyses

- `00.1_runMetaphlan.sh`, `00.2_runHumann3.sh`: run MetaPhlAn & Humann with the latest compatible versions
- `01. Humann_to_GO.sh`: HUMAnN post-processing pipeline. Convert to GO terms.
- `02. GO_term_refinement.py`: GO term post-processing for HUMAnN outputs. Obtain candidate set of GO terms.
- `03. Ancombc2_go_top1000.R`: Differential abundance analysis of HUMAnN gene families mapped to GO terms. 
- `04. Humann_go_ancombc2_plots.R`: Load ANCOM-NC2 results from previous script and generate plots.