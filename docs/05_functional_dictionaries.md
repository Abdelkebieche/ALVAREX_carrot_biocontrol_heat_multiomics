# Functional dictionaries and biological annotation

This repository deliberately separates **statistical recruitment** (WGCNA/maSigPro) from **biological annotation**. A gene first enters a network/trajectory analysis based on expression, then annotation dictionaries are joined at gene level.

## 1. Gene Ontology

GO over-representation is not a keyword dictionary. It uses carrot-specific gene-to-term mappings and analyses BP, MF and CC separately. The **genotype-specific WGCNA retained genes** form the enrichment universe. Significant terms require BH-adjusted P < 0.05 and q < 0.05.

## 2. HSP/chaperone dictionary

The HSP classifier combines explicit annotations and characteristic domains, then excludes known false positives. Families include HSP20/sHSP, HSP40/DnaJ, HSP60/Cpn60, HSP10/Cpn10, CCT/TCP-1, HSP70/HSC70, HSP90, HSP100/ClpB and co-chaperones.

Key safeguards in the supplied implementation include exclusion of GroES-like dehydrogenases, generic HSP70/actin nucleotide-binding-domain matches, generic HATPase kinases/receptors and ambiguous alpha-crystallin catalytic proteins. Subfamilies/localisation labels are assigned only when supported by explicit annotation text; unresolved cases remain conservative. Manual overrides can supersede automated assignments.

## 3. Transcription-factor dictionary

TF identity is annotation-driven rather than inferred from loose function keywords. The companion builder uses genome TF-type/family fields and harmonises interpretation with PlantTFDB/PlantRegMap. DREB is treated as an ERF subgroup when explicitly supported, consistent with AP2/ERF family classification.

## 4. Immune/defence dictionary

The immune classifier is literature-guided and groups genes into receptor-mediated perception, co-receptors/RLCKs, Ca²⁺ signalling, ROS/redox, MAPK signalling, intracellular NLR-associated immunity, defence hormones, transcriptional regulation and downstream defence responses. The classifier records evidence levels and retains a distinction between a broad family assignment and a direct immunity anchor.

A crucial interpretation rule is maintained: heat- or treatment-associated expression of an immune-related gene does **not** by itself demonstrate PTI/ETI activation or altered disease resistance.

## 5. Broad functional families

Eight literature/KEGG-informed families are used:

1. protein folding, processing and degradation;
2. hormone signalling and metabolism;
3. MAPK, calcium and stress signalling;
4. photosynthesis and energy metabolism;
5. terpenoid and specialised metabolism;
6. ROS, redox and glutathione metabolism;
7. lipid and membrane metabolism;
8. carbon, starch, sucrose and osmoprotection.

Membership is built from KEGG provenance, KO/pathway assignments, explicit functional annotations and refined subfamily rules. Supplied subfamily assignments are secondary evidence when stronger curated rules are available. Ambiguous cases are exported for manual review rather than forced into an “Other” biological interpretation.

## Reproducibility dependency

The final biological-heatmap master script reuses annotation outputs created by companion annotation builders. The key builders are included in `scripts/annotation_builders/` as reference implementations; external DH13M14/PlanT2T annotation tables must still be supplied locally.

## Core references

See `references/references.bib` for Liu et al. (2024), Langfelder & Horvath (2008), Yu et al. (2012), Wu et al. (2021), PlantTFDB/PlantRegMap references, and the carrot–*A. dauci* literature used to interpret defence/specialised metabolism.
