#!/bin/bash

set -euo pipefail
command -v STAR >/dev/null || { echo "ERROR: STAR not found in PATH" >&2; exit 1; }

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT before running}"
GENOME_DIR="${GENOME_DIR:-$PROJECT_ROOT/metadata/reference/DH13M14}"
GTF="${GTF_FILE:-$GENOME_DIR/DH13M14.gtf}"
STAR_DIR="${STAR_DIR:-$PROJECT_ROOT/results/03_alignment_or_salmon/STAR}"
OUTDIR="${OUTDIR:-$PROJECT_ROOT/results/04_counts/strand_check_star}"
mkdir -p "$OUTDIR"

BAMS=(
  "$STAR_DIR/0001/0001_Aligned.sortedByCoord.out.bam"
  "$STAR_DIR/0206/0206_Aligned.sortedByCoord.out.bam"
  "$STAR_DIR/0105/0105_Aligned.sortedByCoord.out.bam"
  "$STAR_DIR/0213/0213_Aligned.sortedByCoord.out.bam"
  "$STAR_DIR/0300/0300_Aligned.sortedByCoord.out.bam"
  "$STAR_DIR/0408/0408_Aligned.sortedByCoord.out.bam"
)

for BAM in "${BAMS[@]}"; do
  [[ -f "$BAM" ]] || { echo "ERROR: BAM not found: $BAM" >&2; exit 1; }
  s=$(basename "$(dirname "$BAM")")
  STAR --runThreadN "${THREADS:-10}" --genomeDir "$GENOME_DIR" \
    --sjdbGTFfile "$GTF" --readFilesType SAM PE --inputAlignmentsFromBAM "$BAM" \
    --outFileNamePrefix "$OUTDIR/${s}_" --quantMode GeneCounts --outSAMtype None
done

python3 - "$OUTDIR" <<'PY_STRAND'
import glob, os, sys, csv
outdir=sys.argv[1]
rows=[]
for f in glob.glob(os.path.join(outdir, '*ReadsPerGene.out.tab')):
    sample=os.path.basename(f).split('_')[0]
    sums=[0,0,0]
    with open(f) as h:
        for line in h:
            parts=line.rstrip('\n').split('\t')
            if not parts or parts[0].startswith('N_'):
                continue
            vals=list(map(int, parts[1:4]))
            sums=[a+b for a,b in zip(sums,vals)]
    rows.append((sample,*sums))
rows.sort()
path=os.path.join(outdir,'strand_star_summary.tsv')
with open(path,'w',newline='') as h:
    w=csv.writer(h,delimiter='\t')
    w.writerow(['sample','unstranded_sum','strand1_sum','strand2_sum'])
    w.writerows(rows)
for r in rows: print('\t'.join(map(str,r)))
print('Wrote:', path)
PY_STRAND
