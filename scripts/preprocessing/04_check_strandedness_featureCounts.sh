#!/bin/bash

set -euo pipefail
command -v featureCounts >/dev/null || { echo "ERROR: featureCounts not found in PATH" >&2; exit 1; }

PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT before running}"
STAR_DIR="${STAR_DIR:-$PROJECT_ROOT/results/03_alignment_or_salmon/STAR}"
OUTDIR="${OUTDIR:-$PROJECT_ROOT/results/04_counts/featureCounts/strand_test}"
GTF="${GTF_FILE:-$PROJECT_ROOT/metadata/reference/DH13M14/DH13M14.gtf}"
mkdir -p "$OUTDIR"; cd "$OUTDIR"

BAMS=(
  "$STAR_DIR/0001/0001_Aligned.sortedByCoord.out.bam"
  "$STAR_DIR/0206/0206_Aligned.sortedByCoord.out.bam"
  "$STAR_DIR/0105/0105_Aligned.sortedByCoord.out.bam"
  "$STAR_DIR/0213/0213_Aligned.sortedByCoord.out.bam"
  "$STAR_DIR/0300/0300_Aligned.sortedByCoord.out.bam"
  "$STAR_DIR/0408/0408_Aligned.sortedByCoord.out.bam"
)
[[ -f "$GTF" ]] || { echo "ERROR: GTF not found: $GTF" >&2; exit 1; }
for b in "${BAMS[@]}"; do
  [[ -f "$b" ]] || { echo "ERROR: BAM not found: $b" >&2; exit 1; }
done

for s in 0 1 2; do
  featureCounts -T "${THREADS:-20}" -p -t exon -g gene_id -s "$s" \
    -a "$GTF" -o "test_s${s}.txt" "${BAMS[@]}" > "test_s${s}.log" 2>&1
done

python3 - <<'PY_FC'
import csv, re
from pathlib import Path
out=[]
for s in [0,1,2]:
    p=Path(f'test_s{s}.txt.summary')
    with p.open() as h:
        rows=list(csv.reader(h, delimiter='\t'))
    headers=rows[0][1:]
    data=rows[1:]
    for j,col in enumerate(headers, start=1):
        total=sum(int(r[j]) for r in data)
        assigned=next(int(r[j]) for r in data if r[0]=='Assigned')
        m=re.search(r'/(\d{4})/\1_Aligned\.sortedByCoord\.out\.bam$', col)
        sample=m.group(1) if m else Path(col).name.split('_')[0]
        out.append((s,sample,assigned,total,assigned/total if total else 0))
out.sort(key=lambda x:(x[1],x[0]))
with open('strandness_assigned_rate.tsv','w',newline='') as h:
    w=csv.writer(h,delimiter='\t')
    w.writerow(['strand_s','sample','assigned','total','assigned_rate'])
    w.writerows(out)
for r in out: print('\t'.join(map(str,r)))
print('Wrote: strandness_assigned_rate.tsv')
PY_FC
