#!/bin/bash
#SBATCH --job-name="MetaPhlAn"

module load miniconda/3.10
perl mpa3.pl *1.fastq
