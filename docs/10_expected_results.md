# Expected publication-level checkpoints

These values come from the supplied manuscript and are useful as high-level rerun checks. They are **not** substitutes for exact file checksums.

## Phenotype

- Mean AUDPC Control → Treated: NH 220.27 → 174.31; HS-7 212.00 → 166.21; HS-2 173.42 → 136.85; HS 173.56 → 131.54.
- These correspond to reductions of approximately 20.9%, 21.6%, 21.1% and 24.2%, respectively.
- No significant mean Treatment × Temperature interaction was reported (manuscript P = 0.631).
- Treatment significantly reduced AUDPC in 20 of the 24 Genotype × Temperature cells after FDR correction.
- In untreated plants, the manuscript reports early mean severity of 4.03 (NH), 3.57 (HS-7), 2.71 (HS-2) and 2.83 (HS), and late progression rates of 0.060, 0.121, 0.085 and 0.104 score/day, respectively.

## Transcriptomic WGCNA

- PRESTO: 11,102 genes across 38 D2/D4 samples; 3 non-grey modules.
- ROBILA: 11,317 genes across 37 D2/D4 samples; 8 non-grey modules.
- Soft powers: 12 (PRESTO), 11 (ROBILA).
- Strong treatment association for the brown module in both genotypes.

## maSigPro

- PRESTO: 6,417 selected genes, 8 clusters.
- ROBILA: 1,588 selected genes, 8 clusters.

## Metabolomics

- 17,846 retained UHPLC-HRMS features across 120 samples.
- Reference co-abundance network: 24 ROBILA D4 samples.
- Reference-network soft power: 12.

If a rerun differs substantially, inspect sample exclusions, input versions, annotation mappings, preprocessing choices and package versions before interpreting biological differences.
