# RNA-seq preprocessing and quantification

The manuscript reports strand-specific mRNA PE150 sequencing and the following processing chain:

1. FastQC 0.11.9;
2. fastp 0.23.1 for adapters/low-quality reads and overlapping-read correction;
3. STAR 2.7.9a against the DH13M14 telomere-to-telomere reference;
4. strand-specific featureCounts against the DH13M14 annotation;
5. MultiQC 1.13 for QC aggregation.

Portable shell implementations are in `scripts/preprocessing/`. They contain no scheduler directives or institutional paths.

## Strandedness

The repository retains two diagnostic scripts: STAR `ReadsPerGene.out.tab` and featureCounts `-s 0/1/2`. The final count orientation must be verified from those diagnostics before re-running the publication analysis.

## Important version note

The manuscript reports featureCounts 1.6.1. The portable script invokes `featureCounts` from the user's environment rather than embedding a cluster module version; users reproducing the manuscript should install the reported version or document any deviation.
