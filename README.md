# MHOmicrobiome
 Microbial community analyses on the topic of metabolically healthy obesity. This repository accompanies the source publication [REF].

 

## Table of contents
* [Data](#data)
* [1.- FASTQ file preprocessing](#1-fastq-file-preprocessing)
* [2.- Exploratory metadata analyses](#2-exploratory-metadata-analyses)
* [3.- Exploratory gut microbiome analyses](#3-exploratory-microbiome-analyses)
* [4.- Machine learning scripts](#4-machine-learning-scripts)
* [5.- Network analyses](#5-network-analyses)

## Data
This directory contains supplementary tables 1 and 2:
- `suppTable1.csv`: AI4Food patient metadata & phenotype classification
- `suppTable2.csv`: subject IDs & subject classification obtained from `curatedMetagenomicData`

> Additional patient metadata were obtained from the `curatedMetagenomicData` R package and from study source publications.
> AI4Food raw reads can be accessed at ENA project PRJEBXXXXXX.

## 1 Fastq file preprocessing
Bash & Perl scripts used to preprocess fastq files from the AI4Food project.

## 2 Exploratory metadata analyses
R scripts used to generate tables & supplementary figures related to subject metadata.

## 3 Exploratory gut microbiome analyses
R scripts for batch effect correction, alpha and beta diversity calculation, and differential abundance analyses.

## 4 Machine learning scripts
R scripts used to build random forest, support vector machine and extreme gradient boosting classifiers.

## 5 Network analyses
Scripts in this directory were executed in the following order:
- `create_visualize_networks.R`: co-occurrence network generation & visualization
- `networkMetrics.ipynb`: calculate basic network properties with NetworkX
- `kcores_and_stats.R`: visualize k-cores and calculate statistics for previously computed metrics
- `attacksRandom.ipynb`, `attacksTargeted.ipynb`: perform attacks on networks
