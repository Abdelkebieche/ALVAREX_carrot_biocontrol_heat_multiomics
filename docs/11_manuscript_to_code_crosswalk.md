# Manuscript-to-code crosswalk

This table maps the analysis Methods to repository code. **Original** means the script was present in the supplied analysis archives. **Reconstructed** means the implementation was written from the manuscript specification because the historical script was not supplied. The phenotype implementation is now a supplied original analysis rather than a reconstruction. **Exploratory/legacy** means the code is retained for provenance but is not the final manuscript specification.

| Manuscript section | Analysis | Repository implementation | Status | Key outputs / notes |
|---|---|---|---|---|
| 2.3 | AUDPC, genotype/temperature/treatment mixed models; repeated disease progression | `scripts/phenotype/01_phenotype_analysis.R` | Supplied original analysis; filesystem sanitised | Raw workbook included; full AUDPC/D54 factorial LMM; BH across 24 cell contrasts; AR(1) repeated severity; supplementary early/late block |
| 2.5 | FastQC → fastp → STAR → featureCounts → MultiQC | `scripts/preprocessing/` | Sanitised from supplied workflow | Cluster/scheduler details intentionally removed |
| 2.6 | Genotype-specific RNA WGCNA | `scripts/transcriptomics/wgcna/01–08*` | Original workflow, reorganised | D2+D4; VST; top 50% MAD; signed bicor; PRESTO β12, ROBILA β11 |
| 2.6 | Publication eigengene models | `scripts/transcriptomics/integration/01_compare_PRESTO_ROBILA_publication.R` | Original/final comparative script | `ME ~ Block + Temperature * Treatment * Sampling_time`; Type III; BH |
| 2.7 | GO ORA | `scripts/transcriptomics/wgcna/06_go_enrichment.R` | Original | Genotype-specific network universe |
| 2.7 | Functional dictionaries | `scripts/annotation_builders/` + `docs/05_functional_dictionaries.md` | Reference builders from supplied analysis | HSP, TF, immune/defence, broad functions; ambiguity/manual curation preserved |
| 2.7 | Cross-genotype module comparison | `scripts/transcriptomics/integration/01_compare_PRESTO_ROBILA_publication.R` | Original/final | Directional overlap plus upper-tail hypergeometric tests, BH-adjusted |
| 2.8 | D0/D2/D4 maSigPro | `scripts/transcriptomics/masigpro/01_shared_baseline_maSigPro.R` | Original, portable | Common D0 baseline, degree 2, Q=0.01, 8 clusters |
| 2.8 | WGCNA × maSigPro integration | `scripts/transcriptomics/integration/03_integrate_WGCNA_maSigPro.R` | Original | Hypergeometric enrichment; FDR<0.05, fold enrichment>1, overlap≥10 |
| 2.8 | Biological heatmaps | `scripts/transcriptomics/biological_heatmaps/` | Original workflow, reorganised | D0 shown for biological context; WGCNA membership derives from D2+D4 |
| 2.9 | Compound Discoverer LC-HRMS preprocessing | — | External/vendor workflow not supplied | Raw processing parameters are documented in manuscript but project file was not provided |
| 2.10 | PLS-DA | `scripts/metabolomics/plsda/01_plsda_reconstructed_from_methods.R` | Reconstructed | Per genotype/stage; log2; zero-variance removal; autoscaling; repeated 4-fold CV ×50 |
| 2.10 | LC-MS co-abundance WGCNA | `scripts/metabolomics/wgcna/` | Mixed provenance | Generic D4+D10 workflow is exploratory; manuscript reference network is ROBILA-D4 only |
| 2.11 | Hub-feature projection/manual annotation | `scripts/metabolomics/wgcna/06_project_annotated_hubs_publication.R` | Original final figure logic | Projection ≠ formal module preservation; requires curated feature annotation table |
| supplementary | Temporal metabolite patterns | `scripts/metabolomics/temporal/01_mfuzz_temporal_clustering.R` | Original exploratory/complementary | Mfuzz-based clustering; not the core co-abundance network |

## Why this distinction matters

The repository preserves provenance rather than silently rewriting history. Where the supplied code and manuscript differ, the publication implementation is labelled explicitly and the older script remains available only when it helps audit how the analysis evolved.
