#!/bin/bash   
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10    
#SBATCH --mem=80G             
#SBATCH --time=5:00:00       

# Load R
conda activate qiime2-env
export TMPDIR=./${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$TMPDIR"

# Filenames
FILES=(./subsampling_100iters/*.rds)

# Select filename
FILE=${FILES[$SLURM_ARRAY_TASK_ID-1]}

echo "[$(date)] Processing: $FILE"

# Run R script
Rscript mbScript_subsampledNets.R "$FILE"

rm -rf "$TMPDIR"
