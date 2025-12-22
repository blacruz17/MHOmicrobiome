#!/bin/bash
#SBATCH --job-name=mpa_5M
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=logs/metaphlan_5M_%A_%a.out
#SBATCH --error=logs/metaphlan_5M_%A_%a.err

OUTDIR="./metaphlan_5M"
DB_PATH="./metaphlanJan25"
THREADS_PER_SAMPLE=4

mkdir -p "$OUTDIR" logs
export TMPDIR=/home/tmp_files/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$TMPDIR"

# === SAMPLE SELECTION BASED ON SLURM_ARRAY_TASK_ID ===
# List of bowtie2.bz2 files
BOWTIE_FILES=($(find metaphlan_20250828/bowtie2out -maxdepth 1 -type f -name "*.bowtie2.bz2" | sort))
BOWTIE_FILE="${BOWTIE_FILES[$SLURM_ARRAY_TASK_ID]}"
BASENAME=$(basename "$BOWTIE_FILE" .bowtie2.bz2)
OUT_FILE="${OUTDIR}/${BASENAME}_profile.txt"

# === METAPHLAN ===
echo "[$BASENAME] Processing using $THREADS_PER_SAMPLE threads"
metaphlan "$BOWTIE_FILE" \
    --input_type mapout \
    --output_file "$OUT_FILE" \
    --nproc "$THREADS_PER_SAMPLE" \
    --db_dir "$DB_PATH" \
    --index mpa_vJan25_CHOCOPhlAnSGB_202503 \
    --subsampling 5000000 \
    --mapping_subsampling \
    --subsampling_seed 505 \
    --tmp_dir "$TMPDIR"

rm -rf "$TMPDIR"