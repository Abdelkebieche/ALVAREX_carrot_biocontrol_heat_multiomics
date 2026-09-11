# Bibliography

`references.bib` is the repository-level bibliography for the analysis workflow and biological interpretation.

It is intentionally broader than a software-only citation list. It contains four groups of references:

1. **primary data processing and statistics** — STAR, fastp, featureCounts, MultiQC, lme4/lmerTest, `nlme`, `car`, `emmeans`, BH-FDR;
2. **omics methods** — DESeq2, WGCNA, maSigPro, clusterProfiler, mixOmics and Mfuzz;
3. **annotation/dictionary foundations** — Gene Ontology, KEGG, PlantTFDB/PlantRegMap, AP2/ERF classification, HSP/chaperone family literature and plant-immunity reviews;
4. **carrot–Alternaria / biocontrol context** — DH13M14, carrot resistance QTL, specialised-metabolite studies and Bacillus/QST 2808 studies.

## Bibliographic audit note

The manuscript draft supplied with this repository cites **“Conrath et al. (2006)”** for maSigPro. The canonical maSigPro paper is **Conesa et al. (2006), Bioinformatics 22:1096–1102, DOI 10.1093/bioinformatics/btl056**. The repository uses the canonical citation and flags this as a manuscript citation to correct before submission.

Package versions used for the final analysis should still be captured with `sessionInfo()`/`packageVersion()` because package citations and versions can change independently of the methodological papers.
