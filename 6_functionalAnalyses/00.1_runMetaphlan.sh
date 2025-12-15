FILES_LIST="./fastqFiles.txt" # text file with filenames for all Fastq Files
OUTDIR="./data"
DB_PATH="/home/mpa_vJun23" # path to metaphlan database
THREADS_PER_SAMPLE=4

mkdir -p "$OUTDIR" "$OUTDIR/bowtie2out" logs
export TMPDIR=/home/tmp_files_mpa/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} # path for tmp files
mkdir -p "$TMPDIR"


# === Sample selection (Forward read) based on slurm array ===
FWD=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$FILES_LIST")
# Obtain reverse fastq filename
REV="${FWD/_R1/_R2}"
# Obtain bsaename
BASENAME=$(basename "$FWD" _R1.fq.gz)
OUT_FILE="${OUTDIR}/${BASENAME}_profile.txt"

# === RUN METAPHLAN ===
if [[ -f "$REV" ]]; then
    echo "[$BASENAME] Paired-end using $SLURM_CPUS_PER_TASK threads"
    metaphlan -1 "$FWD" -2 "$REV" \
        --output_file "$OUT_FILE" \
        --subsampling_paired 30000000 \
        --input_type fastq \
        --nproc $SLURM_CPUS_PER_TASK \
	--db_dir "$DB_PATH" \
        --index mpa_vJun23_CHOCOPhlAnSGB_202307 \
        --mapout "${OUTDIR}/bowtie2out/${BASENAME}.bowtie2.bz2" \
	-t rel_ab_w_read_stats \
        --tmp_dir "$TMPDIR"
else
    echo "[$BASENAME] Single-end using $SLURM_CPUS_PER_TASK threads"
    metaphlan "$FWD" \
        --output_file "$OUT_FILE" \
        --input_type fastq \
        --nproc $SLURM_CPUS_PER_TASK \
	--db_dir "$DB_PATH" \
        --index mpa_vJun23_CHOCOPhlAnSGB_202307 \
        --mapout "${OUTDIR}/bowtie2out/${BASENAME}.bowtie2.bz2" \
	-t rel_ab_w_read_stats \
        --tmp_dir "$TMPDIR"
fi

rm -rf "$TMPDIR"
