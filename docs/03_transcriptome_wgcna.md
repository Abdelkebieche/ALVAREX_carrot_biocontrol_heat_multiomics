# Transcriptomic WGCNA

## Publication design

Networks are constructed **separately for PRESTO and ROBILA** from D2 and D4 samples only. Count filtering requires at least 10 total counts and at least 10 counts in three or more samples. Counts are variance-stabilised with DESeq2, zero-variance genes are removed, and the 50% most variable genes by MAD are retained.

WGCNA uses:

- biweight midcorrelation (`bicor`);
- signed adjacency and signed TOM;
- target signed scale-free fit ≈ 0.75;
- publication powers β = 12 (PRESTO) and β = 11 (ROBILA);
- `deepSplit = 3`;
- minimum module size = 100;
- merge cut height = 0.25.

## Module models

The publication model is:

```text
ME ~ Block + Temperature * Treatment * Sampling_time
```

`Block` is a fixed blocking factor. Type III tests use sum-to-zero contrasts. FDR is controlled across modules for each model term with Benjamini-Hochberg. Estimated marginal means are used for within-temperature Treated–Control contrasts and other planned contrasts.

## Hub genes

- initial screen: kME ≥ 0.70;
- high-confidence hub genes: kME ≥ 0.85.

## Script status

`scripts/transcriptomics/wgcna/04_module_trait_exploratory.R` is retained because it was part of the supplied analysis history, but it uses an older simplified model. The publication model is implemented in `scripts/transcriptomics/integration/01_compare_PRESTO_ROBILA_publication.R`.

## Cross-genotype comparison

Only genes retained in both genotype-specific networks are compared. Module correspondence is summarised by shared membership and directional fractions. Hypergeometric overlap tests are corrected by BH. The manuscript's directional figure uses the ROBILA module size as denominator when reporting the fraction mapping to PRESTO modules.
