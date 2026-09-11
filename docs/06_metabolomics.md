# Metabolomics workflow

## Dataset

Untargeted UHPLC-HRMS was analysed for PRESTO and ROBILA at D0, D4 and D10. The manuscript reports 17,846 retained analytical features across 120 samples.

## PLS-DA

The manuscript specifies separate analyses by genotype and stage. Zero-variance features are removed, peak areas are log2-transformed, and variables are autoscaled. At D0 classes are Temperature; at D4/D10 classes are Temperature × Treatment. Performance is evaluated by repeated 4-fold cross-validation with 50 repeats.

The original PLS-DA script was **not included** in the supplied RStudio exports. `scripts/metabolomics/plsda/01_plsda_reconstructed_from_methods.R` is therefore transparently marked as reconstructed from the manuscript and should be checked against the original script before archival release.

## Reference co-abundance network

The publication network is defined on **ROBILA D4 only (24 samples)** so that it shares the D4 stage with RNA-seq while reducing genotype/stage variation. It uses signed bicor, β = 12, `deepSplit = 3`, minimum module size = 100 and merge cut height = 0.25. Module eigengenes are tested with Temperature × Treatment factorial models; BH correction is applied across modules. Grey features are excluded from biological interpretation. Hub features use kME ≥ 0.70.

## Projection, not preservation

Manually annotated hub-feature sets from the ROBILA-D4 reference modules are projected onto PRESTO and ROBILA at D0/D4/D10. The projected score is the mean feature-wise Z-score for the reference hub set. This describes relative abundance of the **same feature set** outside the reference dataset; it is not a formal module-preservation test.

## Important script-history note

The supplied export also contains an exploratory generic LC-MS WGCNA pipeline using D4+D10 and Pareto/MAD filtering. It is retained and explicitly labelled `exploratory`. The publication manuscript instead describes the ROBILA-D4 reference network above. The repository does not silently equate these two analyses.

## Manual feature annotation

Selected hub features were manually examined using accurate monoisotopic mass, isotope pattern, adducts/in-source fragments, MS/MS when available and published specialised-metabolism information. Automated Compound Discoverer labels are not treated as definitive identifications.
