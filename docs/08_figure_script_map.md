# Manuscript figure ↔ script map

| Figure | Analysis | Primary scripts |
|---|---|---|
| V.3 | Disease phenotype | `scripts/phenotype/01_phenotype_analysis.R` (supplied original analysis; public path handling sanitised) |
| V.4 | RNA PCA, WGCNA module associations/eigengenes | `transcriptomics/wgcna/02_*`, `integration/01_compare_PRESTO_ROBILA_publication.R` |
| V.5 | PRESTO ↔ ROBILA module correspondence | `integration/01_compare_PRESTO_ROBILA_publication.R` |
| V.6 | Module functional enrichment summary | `wgcna/06_go_enrichment.R`, `integration/02_functional_summary_publication.R` |
| V.7–V.8 | Defence-associated heatmaps | `biological_heatmaps/00_*`, `04_immune_defence_heatmaps.R` |
| V.9–V.10 | HSP/chaperone heatmaps | `biological_heatmaps/02_hsp_heatmaps.R` |
| V.11 | maSigPro profiles | `masigpro/01_shared_baseline_maSigPro.R` |
| V.12 | WGCNA × maSigPro integration | `integration/03_integrate_WGCNA_maSigPro.R` |
| V.13 | LC-HRMS PLS-DA | reconstructed PLS-DA script; original absent from export |
| V.14 | ROBILA-D4 LC-MS WGCNA | `metabolomics/wgcna/` (publication target documented; supplied generic core retained) |
| V.15–V.16 | Annotated hub heatmaps/projections | `metabolomics/wgcna/06_project_annotated_hubs_publication.R` |
| V.17 | Conceptual synthesis | manuscript schematic; no computational script required |
