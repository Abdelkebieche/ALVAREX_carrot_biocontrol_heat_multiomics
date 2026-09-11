# maSigPro analysis and WGCNA integration

## Shared-baseline design

D0 is pre-treatment, therefore there is **one observed baseline** for each temperature regime. The analysis represents this baseline as shared by the Control and Treated series:

- Control: D0 → D2 → D4
- Treated: D0 (same observed baseline) → D2 → D4

The baseline is duplicated only in the maSigPro design matrix; it is not a second biological observation.

## Filtering and modelling

The supplied script performs count filtering, DESeq2 VST/logCPM generation, then retains the upper 50% of genes by MAD. maSigPro is run separately for each genotype with:

- polynomial degree = 2;
- `p.vector`, BH `Q = 0.01`;
- backward `T.fit`;
- eight hierarchical expression-profile clusters (`ward.D2`).

Because only D0, D2 and D4 are sampled, the clusters should be interpreted as **patterns across observed stages**, not as a dense continuous time course.

## Integration with WGCNA

The publication integration uses strict T.fit genes by default, restricts genes to the corresponding WGCNA universe, and tests maSigPro-cluster × WGCNA-module overlap with an upper-tail hypergeometric test. Highlighted associations require:

- FDR < 0.05;
- fold enrichment > 1;
- ≥ 10 shared genes.

This integration combines two complementary structures: co-expression across D2/D4 and stage-dependent profiles across D0/D2/D4.
