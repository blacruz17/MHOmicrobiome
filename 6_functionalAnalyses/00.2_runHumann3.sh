# 1. List available FASTQ files
FASTQ_FILES=($(ls *.fastq.gz | sort))

# 2. Link basename to FASTQ profile
BASE=$(basename $FASTQ)
SAMPLE=${BASE%.fastq}

# 3. Route to taxonomic profile
TAX=./data/${BASE}_R1.fastq.gz_profile.txt

# 4. Individual output per sample
OUTDIR=HUMANN_OUT/${SAMPLE}
mkdir -p $OUTDIR

# 5. Run HUMAnN
humann \
  -i "$FASTQ" \
  -o "$OUTDIR" \
  --threads 10 \
  --verbose \
  --taxonomic-profile "$TAX" \
  --memory-use maximum
