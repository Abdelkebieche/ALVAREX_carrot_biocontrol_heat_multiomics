#!/bin/bash

set -euo pipefail
for exe in fastqc fastp STAR samtools featureCounts multiqc; do
  command -v "$exe" >/dev/null || { echo "ERROR: $exe not found in PATH" >&2; exit 1; }
done


PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT before running}"
DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/data/raw}"
GENOME_DIR="${GENOME_DIR:-$PROJECT_ROOT/metadata/reference/DH13M14}"
GTF_FILE="${GTF_FILE:-$GENOME_DIR/DH13M14.gtf}"

QC_DIR="$PROJECT_ROOT/results/01_fastqc_multiqc"
TRIM_DIR="$PROJECT_ROOT/results/02_trimming"
MAP_DIR="$PROJECT_ROOT/results/03_alignment_or_salmon"
COUNT_DIR="$PROJECT_ROOT/results/04_counts"
LOG_DIR="$PROJECT_ROOT/results/LOGS"
FASTQC_RAW="$QC_DIR/FastQC_Raw"
FASTQC_CLEAN="$QC_DIR/FastQC_Cleaned"
MULTIQC_RAW="$QC_DIR/MultiQC_Raw"
MULTIQC_CLEAN="$QC_DIR/MultiQC_Cleaned"
TRIM_OUT="$TRIM_DIR/fastp"
STAR_OUT="$MAP_DIR/STAR"
FC_OUT="$COUNT_DIR/featureCounts"

mkdir -p "$FASTQC_RAW" "$FASTQC_CLEAN" "$MULTIQC_RAW" "$MULTIQC_CLEAN" \
  "$TRIM_OUT" "$STAR_OUT" "$FC_OUT" "$LOG_DIR"

[[ -d "$DATA_DIR" ]] || { echo "ERROR: DATA_DIR not found: $DATA_DIR" >&2; exit 1; }
[[ -d "$GENOME_DIR" ]] || { echo "ERROR: GENOME_DIR not found: $GENOME_DIR" >&2; exit 1; }
[[ -f "$GTF_FILE" ]] || { echo "ERROR: GTF not found: $GTF_FILE" >&2; exit 1; }
[[ -f "$GENOME_DIR/SA" ]] || echo "WARNING: STAR index marker $GENOME_DIR/SA not found" >&2

CPUS="${THREADS:-20}"
FASTQC_T=$(( CPUS < 20 ? CPUS : 20 ))
FASTP_T=$(( CPUS < 20 ? CPUS : 20 ))
STAR_T="$CPUS"
SAM_T=$(( CPUS < 8 ? CPUS : 8 ))
FC_T=$(( CPUS < 20 ? CPUS : 20 ))

echo "Threads: FastQC=$FASTQC_T fastp=$FASTP_T STAR=$STAR_T samtools=$SAM_T featureCounts=$FC_T"

shopt -s nullglob
SAMPLE_DIRS=( "$DATA_DIR"/*/ )
(( ${#SAMPLE_DIRS[@]} > 0 )) || { echo "ERROR: no sample directories in $DATA_DIR" >&2; exit 1; }

for sample_dir in "${SAMPLE_DIRS[@]}"; do
  sample="$(basename "$sample_dir")"
  BAM_DONE="$STAR_OUT/$sample/${sample}_Aligned.sortedByCoord.out.bam"
  if [[ -f "$BAM_DONE" ]]; then
    echo "Already done: $sample -> skip"
    continue
  fi

  R1_raw="$sample_dir/${sample}_1.fq.gz"
  R2_raw="$sample_dir/${sample}_2.fq.gz"
  [[ -f "$R1_raw" ]] || R1_raw="$sample_dir/${sample}_1.fastq.gz"
  [[ -f "$R2_raw" ]] || R2_raw="$sample_dir/${sample}_2.fastq.gz"
  [[ -f "$R1_raw" ]] || R1_raw="$sample_dir/${sample}_R1.fastq.gz"
  [[ -f "$R2_raw" ]] || R2_raw="$sample_dir/${sample}_R2.fastq.gz"
  [[ -f "$R1_raw" ]] || R1_raw="$sample_dir/${sample}_R1.fq.gz"
  [[ -f "$R2_raw" ]] || R2_raw="$sample_dir/${sample}_R2.fq.gz"
  if [[ ! -f "$R1_raw" || ! -f "$R2_raw" ]]; then
    echo "WARNING: paired reads missing for $sample -> skip" >&2
    continue
  fi

  mkdir -p "$FASTQC_RAW/$sample" "$FASTQC_CLEAN/$sample" "$TRIM_OUT/$sample" \
    "$STAR_OUT/$sample" "$LOG_DIR/$sample"

  fastqc -t "$FASTQC_T" -q --outdir "$FASTQC_RAW/$sample" "$R1_raw" "$R2_raw" \
    > "$LOG_DIR/$sample/fastqc_raw.log" 2>&1

  R1_clean="$TRIM_OUT/$sample/${sample}_clean_1.fq.gz"
  R2_clean="$TRIM_OUT/$sample/${sample}_clean_2.fq.gz"
  fastp -i "$R1_raw" -I "$R2_raw" -o "$R1_clean" -O "$R2_clean" \
    -w "$FASTP_T" --html "$TRIM_OUT/$sample/${sample}_fastp.html" \
    --json "$TRIM_OUT/$sample/${sample}_fastp.json" \
    > "$LOG_DIR/$sample/fastp.log" 2>&1

  fastqc -t "$FASTQC_T" -q --outdir "$FASTQC_CLEAN/$sample" "$R1_clean" "$R2_clean" \
    > "$LOG_DIR/$sample/fastqc_clean.log" 2>&1

  STAR_PREFIX="$STAR_OUT/$sample/${sample}_"
  STAR --genomeDir "$GENOME_DIR" --runThreadN "$STAR_T" --readFilesCommand zcat \
    --readFilesIn "$R1_clean" "$R2_clean" --outFileNamePrefix "$STAR_PREFIX" \
    --outSAMtype BAM SortedByCoordinate --sjdbGTFfile "$GTF_FILE" \
    --alignIntronMax 60000 --alignMatesGapMax 10000 --outFilterMultimapNmax 10 \
    --limitBAMsortRAM 20000000000 --outSAMstrandField intronMotif \
    > "$LOG_DIR/$sample/STAR.log" 2>&1

  BAM="${STAR_PREFIX}Aligned.sortedByCoord.out.bam"
  [[ -f "$BAM" ]] || { echo "ERROR: STAR BAM missing for $sample: $BAM" >&2; exit 1; }
  samtools index -@ "$SAM_T" "$BAM"
  samtools flagstat -@ "$SAM_T" "$BAM" > "$LOG_DIR/$sample/flagstat.txt"
done

multiqc "$FASTQC_RAW" --dirs -o "$MULTIQC_RAW" > "$LOG_DIR/multiqc_raw.log" 2>&1
multiqc "$FASTQC_CLEAN" --dirs -o "$MULTIQC_CLEAN" > "$LOG_DIR/multiqc_clean.log" 2>&1

BAMS=( "$STAR_OUT"/*/*Aligned.sortedByCoord.out.bam )
(( ${#BAMS[@]} > 0 )) || { echo "ERROR: no BAM files found in $STAR_OUT" >&2; exit 1; }

# Final counting setting retained from the original analysis: reverse-stranded (-s 2).
featureCounts -T "$FC_T" -p -t exon -g gene_id -s 2 \
  -a "$GTF_FILE" -o "$FC_OUT/gene_counts_featureCounts.txt" "${BAMS[@]}" \
  > "$LOG_DIR/featureCounts.log" 2>&1

echo "Pipeline completed"
echo "Counts: $FC_OUT/gene_counts_featureCounts.txt"
echo "Logs:   $LOG_DIR"
