# Functional dictionary reference builders

These scripts preserve the detailed keyword/domain/annotation rules used to curate the HSP/chaperone, transcription-factor, immune/defence and broad functional dictionaries.

They are included primarily for **method transparency and auditability**. They originate from the preceding heat-response analysis and are reused here as annotation logic; they are not the statistical engine of the WGCNA article.

## Portable compatibility layer

Run the builders from this directory. `00_builder_config.R` maps the historical variable names onto the repository-wide portable configuration in `../../config/config.R`; `00_builder_helpers.R` provides the shared helper functions. No HPC or institutional paths are stored.

The builders may require external annotation resources and, for figures that reproduce the earlier differential-expression curation, intermediate edgeR result tables. Those resources are deliberately not fabricated in this repository. See `../../docs/05_functional_dictionaries.md` and `../../data/annotation/README.md`.

## Scope

- `01_hsp_dictionary_builder.R`: conservative HSP/chaperone family assignment and exclusion rules.
- `02_tf_dictionary_builder.R`: transcription-factor family harmonisation and ambiguity checks.
- `03_immune_dictionary_builder.R`: literature-guided immune/defence repertoire with confidence tiers.
- `04_functional_dictionary_curator.R`: broad biological programmes and curated subfamilies.

The final WGCNA article uses these curated biological categories for module interpretation; network statistics are calculated independently of dictionary membership.
