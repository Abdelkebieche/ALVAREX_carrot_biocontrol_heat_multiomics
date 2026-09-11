#!/usr/bin/env python3
"""Static repository checks that do not require R or the raw omics data."""
from pathlib import Path
import csv, re, sys

ROOT = Path(__file__).resolve().parents[2]
errors = []

# 1. Security scan: do not publish cluster/institutional execution details.
patterns = {
    "absolute shared path": r"/shared/",
    "scheduler directive": r"#SBATCH",
    "SLURM variable": r"SLURM_",
    "institutional notification email": r"@inrae\.fr",
}
text_ext = {'.R','.r','.sh','.md','.txt','.tsv','.csv','.yml','.yaml','.cff','.bib','.py'}
for p in ROOT.rglob('*'):
    if not p.is_file() or p.suffix.lower() not in text_ext or p.name in {'MANIFEST.sha256.tsv', Path(__file__).name}:
        continue
    txt = p.read_text(errors='ignore')
    for label, pattern in patterns.items():
        if re.search(pattern, txt, flags=re.I|re.M):
            errors.append(f"{label}: {p.relative_to(ROOT)}")

# 2. Literal source() targets: scripts are run from their containing directory.
source_re = re.compile(r"source\(\s*['\"]([^'\"]+)['\"]\s*\)")
for p in ROOT.rglob('*.R'):
    txt = p.read_text(errors='ignore')
    for src in source_re.findall(txt):
        target = (p.parent / src).resolve()
        if not target.exists():
            errors.append(f"missing source target: {p.relative_to(ROOT)} -> {src}")

# 3. RNA metadata integrity for the supplied transcriptomic layer.
meta = ROOT / 'metadata' / 'rnaseq_metadata.csv'
if meta.exists():
    with meta.open(newline='', encoding='utf-8-sig') as fh:
        sample = fh.read(4096)
        delim = ';' if sample.count(';') > sample.count(',') else ','
        fh.seek(0)
        rows = list(csv.DictReader(fh, delimiter=delim))
    required = {'Sample','Genotype','Temperature','Treatment','Block','Sampling_time'}
    cols = set(rows[0]) if rows else set()
    missing = required - cols
    if missing:
        errors.append('metadata missing columns: ' + ', '.join(sorted(missing)))
    genotypes = {r.get('Genotype','').strip() for r in rows}
    if not {'PRESTO','ROBILA'}.issubset(genotypes):
        errors.append('RNA metadata does not contain both PRESTO and ROBILA')
else:
    errors.append('missing metadata/rnaseq_metadata.csv')

if errors:
    print('REPOSITORY VALIDATION FAILED')
    for e in errors:
        print(' -', e)
    sys.exit(1)

print('REPOSITORY VALIDATION PASSED')
print(' - no cluster/HPC disclosure patterns detected')
print(' - all literal R source() targets exist')
print(' - RNA metadata schema is consistent with PRESTO/ROBILA omics layer')
