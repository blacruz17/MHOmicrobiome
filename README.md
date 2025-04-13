# MHOmicrobiome
Microbial community analyses on the topic of metabolically healthy obesity.

This repository contains a set of scripts used to preprocess sequencing data, analyze metadata, explore gut microbiome profiles, build machine learning models, and perform network-based analyses. Our code focuses on the comparison of four groups encompassing different metabolic health and obesity phenotypes.

## Table of contents
* [1.- FASTQ file preprocessing](#1-fastq-file-preprocessing)
* [2.- Exploratory metadata analyses](#2-exploratory-metadata-analyses)
* [3.- Exploratory gut microbiome analyses](#3-exploratory-gut-microbiome-analyses)
* [4.- Machine learning scripts](#4-machine-learning-scripts)
* [5.- Network analyses](#5-network-analyses)

## 1 Fastq file preprocessing

Bash & Perl scripts used to preprocess fastq files from the AI4Food project. This folder contains the following scripts:

- `metawrap.sh`: performs quality control as well as sequence filtering and trimming. We did this on metaWRAP version 1.2.1 using the hg38 version of the human genome.
- `mpa3.pl`: uses MetaPhlAn for taxonomic assignments. This code uses MetaPhlAn version 4.1 with 
the “mpa3” parameter and maps against the mpa_v30_CHOCOPhlAn_201901 database. This script was executed using the `slurm_mpa3.sh` scritpt.

No preprocessing of the remaining datasets was performed since the taxonomic assignments were already available within the `curatedMetagenomicData` R package.

## 2 Exploratory metadata analyses
R scripts used to generate tables & supplementary figures related to subject metadata:

- `metadataExploration.R`: compares biochemical and anthropometrical metadata among the four groups, generating a table with all comparisons (Table 2).
- `supp_heatmap.R`: explores all datasets within `curatedMetagenomicData` to look for studies with stool samples from adult subjects not taking antibiotics, and where at least one of our variables of interest to describe obesity and metabolic phenotypes was availabl. Generates Supplementary Figure 1.
- `supp_boxplots.R`: in-depth statistical analysis accompanying the previous scripts as well as the necessary code to replicate Supplementary Figure 2.

Library requirements are specified in the beggining of each R script.

## 3 Exploratory gut microbiome analyses
R scripts for batch effect correction, alpha and beta diversity calculation, and differential abundance analyses:

- `mmuphin_batch_correction.R`: this is the first of the 3 scripts to be executed. It performs batch effect correction caused by the study of origin while controlling for the effect of metabolic health and obesity using MMUPHin.
- `diversity.R`: alpha- and beta-diversity exploration. Generates Figures 1A and 1B.
- `diffAbundance.R`: uses Analysis of Compositions of Microbiomes with Bias Correction 2 (ANCOM-BC2), an extension of the ANCOM-BC methodology that can be implemented in datasets with multiple groups. Controls for multiple pairwise comparisons based on the mixed directional false discovery rate (mdFDR) using the Holm-Bonferroni procedure. Generates Figure 1C.

Library requirements are specified in the beggining of each R script.

## 4 Machine learning scripts
R scripts used to build random forest, support vector machine and extreme gradient boosting classifiers.


Library requirements are specified in the beggining of each R script.

## 5 Network analyses
Scripts in this directory were executed in the following order:
- `create_visualize_networks.R`: co-occurrence network generation & visualization
- `networkMetrics.ipynb`: calculate basic network properties with NetworkX
- `kcores_and_stats.R`: visualize k-cores and calculate statistics for previously computed metrics
- `attacksRandom.ipynb`, `attacksTargeted.ipynb`: perform attacks on networks

These scripts can be run locally and can also be easily executed in cloud services such as Google Colab.
