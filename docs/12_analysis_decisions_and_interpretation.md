# Analysis decisions and interpretation guardrails

These rules are part of the reproducible analysis because they determine what the outputs can and cannot support biologically.

## Experimental scope

- Six genotypes are used for **disease phenotyping**.
- PRESTO and ROBILA are used for **RNA-seq and UHPLC-HRMS**.
- NH, HS-7 and HS-2 are used for the omics analyses; the HS condition is retained in the phenotype but excluded from omics because heat started at D0, simultaneously with treatment.
- HS-7 and HS-2 are distinct thermal histories, not successive time points of one heat-duration series.

## Phenotype inference

- AUDPC is the primary disease endpoint; D54 is secondary.
- Cell-wise Control-vs-Treated inference uses 24 Temperature × Genotype contrasts from the full factorial AUDPC model, with BH correction across the 24 tests.
- The repeated-measures AR(1) model is used to characterise longitudinal disease progression rather than to generate significance tests at every assessment date.
- The supplied early/late supplementary block tests temperature effects in Control plants only. The manuscript currently describes a broader Treatment-inclusive model for these derived variables. This difference must be resolved before final submission.
- Descriptive percentage reduction in AUDPC is not itself a test of Treatment × Genotype or Treatment × Temperature interaction.

## RNA-seq WGCNA

- Networks are genotype-specific and use D2+D4 only.
- D0 is not allowed to influence WGCNA module construction; it is added later only for selected visualisations.
- Biological replicate values, not group means, are used for statistical models.
- Group means and gene-wise Z-scores are used for display only.
- Module labels such as `brown`, `blue`, or `turquoise` are arbitrary within each independently constructed network. Cross-genotype biological equivalence is tested from shared genes, not inferred from matching colour names alone.

## maSigPro shared baseline

At D0, treatment has not yet been applied. The Control D0 measurements are therefore used as the common starting state for both the Control and Treated series. Duplicating D0 values in the design is a modelling device; it does **not** create additional biological replicates.

With only D0, D2 and D4, the fitted polynomial profiles describe patterns across the observed stages. They should not be interpreted as densely sampled kinetic trajectories.

## D2 → D4 interpretation

D2 is immediately before *A. dauci* inoculation and D4 is two days after inoculation. Because no matched non-inoculated D4 control was included, a D2-to-D4 change cannot be attributed uniquely to pathogen challenge. It may include stage, continued thermal exposure, treatment dynamics and infection-related effects.

## WGCNA × maSigPro

Cluster–module correspondence is tested with an upper-tail hypergeometric test. Reported associations require all of:

- FDR < 0.05;
- fold enrichment > 1;
- at least 10 shared genes.

This is an enrichment/cross-classification analysis, not a causal link between a maSigPro pattern and a WGCNA module.

## Functional dictionaries

Dictionary membership supports interpretation; it does not determine network construction. Ambiguous annotations remain ambiguous or are manually curated rather than forced into a category. Expression of an immune-associated gene is not equivalent to activation of PTI/ETI, and an HSP/defence label is not evidence of altered resistance by itself.

## Metabolomics

- PLS-DA is treated as supervised exploratory visualisation and must be accompanied by repeated cross-validation.
- The manuscript reference co-abundance network is ROBILA at D4.
- Hub features are defined within that reference network by kME.
- Projecting the same hub-feature set into PRESTO, D0 or D10 measures relative abundance of those features; it is **not** a formal module-preservation test.
- Manual compound annotation should report the evidence available (accurate mass, isotope pattern, adducts, in-source fragments, MS/MS and literature) and should not upgrade tentative identities beyond the evidence.

## Cross-omics interpretation

Transcriptomic and metabolomic programmes are compared at the level of independently identified functions/programmes. No direct gene–metabolite association is claimed unless an explicit statistical test is added in a future analysis.
