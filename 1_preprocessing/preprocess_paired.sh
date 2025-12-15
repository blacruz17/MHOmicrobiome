#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --output=/home/FengQ_2015/paired/log/preprocess.log
#SBATCH --error=/home/FengQ_2015/paired/log/preprocess.log
#SBATCH --cpus-per-task=30
#SBATCH --time=48:00:00

export PREPRO=/home/Nextflow-scripts

nextflow run $PREPRO/preprocess.nf \
	--cpus 30 \
	--data_dir /home/FengQ_2015/paired \
	--human_index /home/Homosapiens_genome/Bowtie2_index \
