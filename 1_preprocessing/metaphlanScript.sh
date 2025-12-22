#!/bin/bash
#SBATCH --job-name=metaphlan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=logs/metaphlan_%A_%a.out
#SBATCH --error=logs/metaphlan_%A_%a.err


FILES_LIST="./fastqFiles.txt"
OUTDIR="./metaphlan_$(date +%Y%m%d)"
DB_PATH="./metaphlanJan25"
THREADS_PER_SAMPLE=4

mkdir -p "$OUTDIR" "$OUTDIR/bowtie2out" logs

export TMPDIR=/tmp_files

# === SAMPLE SELECTION BASED ON SLURM_ARRAY_TASK_ID ===
FWD=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$FILES_LIST")
REV="${FWD/_R1/_R2}"
BASENAME=$(basename "$FWD" _R1.fq.gz)
OUT_FILE="${OUTDIR}/${BASENAME}_profile.txt"


# === METAPHLAN ===
if [[ -f "$REV" ]]; then
    echo "[$BASENAME] Paired-end using $SLURM_CPUS_PER_TASK threads"
    metaphlan -1 "$FWD" -2 "$REV" \
        --subsampling_paired 30000000 \
        --output_file "$OUT_FILE" \
        --input_type fastq \
        --nproc $SLURM_CPUS_PER_TASK \
	--db_dir "$DB_PATH" \
        --index mpa_vJan25_CHOCOPhlAnSGB_202503 \
        --mapout "${OUTDIR}/bowtie2out/${BASENAME}.bowtie2.bz2"
else
    echo "[$BASENAME] Single-end using $SLURM_CPUS_PER_TASK threads"
    metaphlan "$FWD" \
        --output_file "$OUT_FILE" \
        --input_type fastq \
        --nproc $SLURM_CPUS_PER_TASK \
        --db_dir "$DB_TO_USE" \
        --index mpa_vJan25_CHOCOPhlAnSGB_202503 \
        --mapout "${OUTDIR}/bowtie2out/${BASENAME}.bowtie2.bz2"

fi

# clean temp files
rm -rf $TMPDIR