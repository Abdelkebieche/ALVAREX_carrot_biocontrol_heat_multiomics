# Phenotype analysis module

## Source and provenance

`01_phenotype_analysis.R` is the phenotype-analysis script supplied with the project. In this public repository, only filesystem handling was changed: the institutional path was removed, the input is read through `config/config.R`, and outputs are written below `outputs/phenotype/`. The statistical analysis itself was not rewritten.

## Run

From the repository root:

```bash
Rscript scripts/phenotype/01_phenotype_analysis.R
```

The default input is `data/phenotype/Data_phenotype.xlsx`. Override it without editing tracked code with:

```bash
export CARROT_PHENOTYPE_FILE=/path/to/Data_phenotype.xlsx
export CARROT_OUTPUT_ROOT=/path/to/outputs
Rscript scripts/phenotype/01_phenotype_analysis.R
```

## Main statistical workflow

1. Recode Water/PRI as Control/Treated and T0/T1/T2/T3 as NH/HS-7/HS-2/HS.
2. Compute trapezoidal AUDPC from the seven disease assessments (11–54 DAI).
3. Fit the primary mixed model:

   `AUDPC ~ Treatment * Temperature * Genotype + (1 | Rep)`

4. Fit the secondary D54 model with the same factorial structure.
5. Obtain Type-III tests using sum-to-zero contrasts.
6. Compute Control-vs-Treated `emmeans` contrasts within each Temperature × Genotype cell and apply Benjamini–Hochberg correction across the 24 cells.
7. Compute Tukey compact-letter displays for Temperature and Genotype marginal means.
8. Fit an `nlme::lme` repeated-measures model with AR(1) residual correlation for descriptive disease progression.
9. Check AUDPC/D54 residuals using Shapiro–Wilk, Anderson–Darling and DHARMa diagnostics.
10. Analyse early (11/18/26 DAI) and late (32/40/47/54 DAI) temperature-associated disease dynamics in the supplementary two-phase block.

## Important manuscript-alignment note

The supplied two-phase block performs the **formal early-level and late-rate temperature tests on Control plants only** (`Temperature * Genotype + (1 | Rep)`). The current manuscript Methods state that early level and late progression were analysed with Treatment, Temperature, Genotype and their interactions. These two descriptions are not identical. Before final submission, either the Methods or the analysis implementation should be aligned and the final choice documented. The repository does not silently alter this distinction.

## Outputs

The script writes to `outputs/phenotype/`:

- `figures/pdf`, `figures/png`, `figures/svg`;
- `tables/T1_...T11_...`;
- `tables/Ttp1_...Ttp7_...` for the early/late supplementary block;
- `diagnostics/DHARMa_AUDPC.pdf`.

See `docs/07_phenotype.md` for the full method and references.
