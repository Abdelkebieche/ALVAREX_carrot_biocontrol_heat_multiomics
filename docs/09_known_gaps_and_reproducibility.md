# Reproducibility audit and known gaps

This file records differences between the manuscript description and the supplied script archive. These are intentionally visible.

## 1. Phenotype source code recovered

The supplied phenotype-analysis R script and its original phenotype workbook are now included. The earlier Methods-based reconstruction has been removed. Public-repository changes are limited to portable filesystem/output handling.

### Phenotype manuscript/code alignment still required

The script's supplementary early/late formal models are restricted to Control plants (`Temperature * Genotype + (1|Rep)`), whereas the current manuscript Methods state that Treatment, Temperature, Genotype and their interactions were included for early and late derived variables. This must be reconciled before final publication. The repository preserves the supplied implementation and makes the discrepancy explicit.

## 2. PLS-DA source code absent

The supplied LC-MS export did not contain a `mixOmics` PLS-DA script. A reconstructed implementation based on the manuscript specification is included and labelled as reconstructed.

## 3. Transcriptomic module model has exploratory and publication variants

An older WGCNA script fits `ME ~ Temperature * Treatment + Time` with Type II tests. The final comparative script fits the manuscript model `ME ~ Block + Temperature * Treatment * Sampling_time` with Type III tests and sum contrasts. The latter is the publication model.

## 4. LC-MS network scope differs between generic and publication descriptions

The generic supplied LC-MS workflow filters P3+P4 (D4+D10) and uses log2/Pareto/MAD preprocessing. The manuscript describes a reference co-abundance network built on ROBILA D4 only. The final annotated-hub script points to an historical ROBILA-P3 WGCNA run not included in the archive. Before final public release, the exact script/configuration that created that reference run should be added if available.

## 5. Annotation resources not bundled

DH13M14 functional annotation, PlanT2T gene-ID mapping and some intermediate dictionary outputs are external dependencies. Their expected locations are generic and local-only in `config/config.R`.

## 6. No cluster information is published

All institutional HPC paths, allocations, scheduler directives, partitions, notification emails and private project-directory names have been removed.

## 7. maSigPro citation in the manuscript draft

The draft cites “Conrath et al. (2006)” for maSigPro. The canonical maSigPro method paper is Conesa et al. (2006; DOI 10.1093/bioinformatics/btl056). The repository bibliography uses the canonical citation; the manuscript citation should be checked before submission.

## 8. featureCounts/Subread version requires confirmation

The manuscript draft reports featureCounts v1.6.1, whereas the supplied historical HPC pipeline loaded Subread v2.0.6. Because this can matter for exact reproducibility, the final public release should confirm which executable produced the count matrix used in the manuscript and report that exact version consistently in the Methods and `environment/software_versions.tsv`.
