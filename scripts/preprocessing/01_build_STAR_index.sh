#!/bin/bash

set -euo pipefail
command -v STAR >/dev/null || { echo "ERROR: STAR not found in PATH" >&2; exit 1; }

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT before running}"
REF_DIR="${REF_DIR:-$PROJECT_ROOT/metadata/reference/DH13M14}"
GENOME_FASTA="${GENOME_FASTA:-$REF_DIR/DH13M14_genome.fasta}"
GTF_FILE="${GTF_FILE:-$REF_DIR/DH13M14.gtf}"
THREADS="${THREADS:-24}"

mkdir -p "$REF_DIR"
[[ -f "$GENOME_FASTA" ]] || { echo "ERROR: genome FASTA not found: $GENOME_FASTA" >&2; exit 1; }
[[ -f "$GTF_FILE" ]] || { echo "ERROR: GTF not found: $GTF_FILE" >&2; exit 1; }

STAR --runThreadN "$THREADS" \
  --runMode genomeGenerate \
  --genomeDir "$REF_DIR" \
  --genomeFastaFiles "$GENOME_FASTA" \
  --sjdbGTFfile "$GTF_FILE" \
  --sjdbOverhang 99 \
  --genomeSAindexNbases 13

echo "STAR index written to: $REF_DIR"
