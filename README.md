# MHOmicrobiome
Microbial community analyses on the topic of metabolically healthy obesity.

This repository contains a set of scripts used to preprocess sequencing data, analyze metadata, explore gut microbiome profiles, build machine learning models, and perform network-based analyses. Our code focuses on the comparison of four groups encompassing different metabolic health and obesity phenotypes.

## Table of contents
* [Data](#data)
* [1.- FASTQ file preprocessing](#1-fastq-file-preprocessing)
* [2.- Exploratory metadata analyses](#2-exploratory-metadata-analyses)
* [3.- Exploratory gut microbiome analyses](#3-exploratory-gut-microbiome-analyses)
* [4.- Machine learning scripts](#4-machine-learning-scripts)
* [5.- Network analyses](#5-network-analyses)


## Data 
We provide the following files for those interested in replicating our analyses:

- `physeqMHO_raw.rds`: phyloseq object before batch correction
- Network files: for each phenotype, a file with edge properties (e.g.: `mho_edges.csv`) and a file with node properties (e.g.: `mho_metadata.csv`) are given.


## 1 Fastq file preprocessing

Bash & Perl scripts used to preprocess fastq files from the AI4Food project. This folder contains the following scripts:

- `metawrap.sh`: performs quality control as well as sequence filtering and trimming. We did this on metaWRAP version 1.2.1 using the hg38 version of the human genome.
- `mpa3.pl`: uses MetaPhlAn for taxonomic assignments. This code uses MetaPhlAn version 4.1 with 
the “mpa3” parameter and maps against the mpa_v30_CHOCOPhlAn_201901 database. This script was executed using the `slurm_mpa3.sh` scritpt.

No preprocessing of the remaining datasets was performed since the taxonomic assignments were already available within the `curatedMetagenomicData` R package (version 3.8).

To replicate these analyses, download the AI4Food reads from ENA accession number PRJEB87701.

## 2 Exploratory metadata analyses
R scripts used to generate tables & supplementary figures related to subject metadata:

- `metadataExploration.R`: compares biochemical and anthropometrical metadata among the four groups, generating a table with all comparisons (Table 2).
- `supp_heatmap.R`: explores all datasets within `curatedMetagenomicData` (version 3.8) to look for studies with stool samples from adult subjects not taking antibiotics, and where at least one of our variables of interest to describe obesity and metabolic phenotypes was available. Generates Supplementary Figure 1.
- `supp_boxplots.R`: in-depth statistical analysis accompanying the previous scripts as well as the necessary code to replicate Supplementary Figure 2.

In order to fully recapitulate our analyses, you may need:
- Supplementary Table 1 with AI4Food patient metadata.
- Supplementary files from `curatedMetagenomicData` datasets' source publications
- Sample metadata already available in `curatedMetagenomicData`

Library requirements are specified in the beggining of each R script.

## 3 Exploratory gut microbiome analyses
R scripts for batch effect correction, alpha and beta diversity calculation, and differential abundance analyses:

- `mmuphin_batch_correction.R`: this is the first of the 3 scripts to be executed. It performs batch effect correction caused by the study of origin while controlling for the effect of metabolic health and obesity using MMUPHin (version 1.14).
- `diversity.R`: alpha- and beta-diversity exploration. Generates Figures 1A and 1B. Uses R libraries microbiome (version 1.22) and phyloseq (version 1.46).
- `diffAbundance.R`: uses Analysis of Compositions of Microbiomes with Bias Correction 2 (ANCOM-BC2), an extension of the ANCOM-BC methodology that can be implemented in datasets with multiple groups (version 2.2.2). Controls for multiple pairwise comparisons based on the mixed directional false discovery rate (mdFDR) using the Holm-Bonferroni procedure. Generates Figure 1C.

Library requirements are specified in the beggining of each R script.

## 4 Machine learning scripts
R scripts used to build random forest, support vector machine and extreme gradient boosting classifiers:

- `randomForest/`: contains scripts `basicModels.R` for random forest models without subsampling and `subsamplingModels.R` for models using upsampling or downsampling. Uses the ranger library (version 0.17) 
- `svm/`: contains scripts for basic SVMs (`svm_basic.R`), SVMs employing weighting and/or subsampling approaches (`svm_weighted_subsampling.R`), and multiclass SVMs (`svm_multiclass.R`). Requires kernlab version 0.9-33.
- `xgboost/`: contains scripts used for hyperparameter tuning ( `xgb_hyperparameters_step1.R`, `xgb_hyperparameters_step2.R`) as well as model training (`xgb_train.R`). Requires xgboost version 1.7.8.1.
- Scripts for results analysis: `binaryROCs.R` (only binary models) `multi_and_binaryROCs.R` (multiclass models and Figure 2 creation).

All models were built using the caret R package (version 6.0-94). Further library requirements are specified in the beggining of each R script.

## 5 Network analyses
Scripts in this directory were executed in the following order:
- `create_visualize_networks.R`: co-occurrence network generation & visualization using NetCoMi (version 1.1) and SpiecEasi (version 1.1.3) R libraries. Generates Figures 3A-D. Exports networks as csv files for their use in subsequent analyses.
- `networkMetrics.ipynb`: calculates basic network properties with NetworkX (version 3.4.2) for Table 3.
- `kcores_and_stats.R`: visualize k-cores, generates Figure 3E, and calculates statistics for previously computed metrics (Table 3)
- `attacksRandom.ipynb`, `attacksTargeted.ipynb`: these scritps perform attacks on networks:
  - random attacks, where 1000 randomizations are performed, using a seed for reproducibility
  - targeted attacks based on degree, betweenness, and mean abundance within the dataset.

These scripts can be run locally and can also be easily executed in cloud services such as Google Colab.
