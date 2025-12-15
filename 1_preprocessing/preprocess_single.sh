#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --output=/home/Metacardis/preprocess.log
#SBATCH --error=/home/Metacardis/preprocess.log
#SBATCH --cpus-per-task=30
#SBATCH --time=148:00:00

export PREPRO=/home/Nextflow-scripts

nextflow run $PREPRO/preprocess.nf \
	--cpus 30 \
	--data_dir /home/Metacardis \
	--human_index /home/Homosapiens_genome/Bowtie2_index \
	--bac_filtering false \
	--single_end true
