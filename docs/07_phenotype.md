# Phenotype analysis

## Scope

The phenotype layer uses all six carrot genotypes (`DEEP`, `NANT`, `NEVA`, `OXHELLA`, `PRESTO`, `ROBILA`) and all four thermal regimes (`NH`, `HS-7`, `HS-2`, `HS`). The supplied workbook contains 192 experimental units, corresponding to 2 Treatments × 4 Temperatures × 6 Genotypes × 4 Replicates/blocks. Disease was assessed at 11, 18, 26, 32, 40, 47 and 54 days after inoculation (DAI).

The original analysis is now included as `scripts/phenotype/01_phenotype_analysis.R`; the earlier repository reconstruction has been removed. The public version changes only path/output handling.

## Input and recoding

The original workbook is `data/phenotype/Data_phenotype.xlsx`. The script recodes:

| Raw value | Analysis label |
|---|---|
| Water | Control |
| PRI | Treated |
| T0 | NH |
| T1 | HS-7 |
| T2 | HS-2 |
| T3 | HS |

The disease response columns are D11, D18, D26, D32, D40, D47 and D54. Biomass columns present in the workbook are not retained in the disease-analysis dataframe.

## Primary response: AUDPC

For each experimental unit, AUDPC is calculated by trapezoidal integration:

`AUDPC = Σ (t[i+1] - t[i]) × (y[i] + y[i+1]) / 2`

across the seven disease assessments. The script also calculates relative AUDPC as a percentage of the theoretical maximum score-time area, but absolute AUDPC is the primary inferential response.

The primary model is:

```text
AUDPC ~ Treatment * Temperature * Genotype + (1 | Rep)
```

The model is fitted with `lme4::lmer` using REML. Sum-to-zero contrasts (`contr.sum`) are set before Type-III tests so that factorial main effects/interactions are interpretable in the usual Type-III sense. Type-III tests are obtained with `lmerTest`.

### Treatment contrasts

Estimated marginal means are calculated using `emmeans`. Control and Treated are contrasted separately within each Temperature × Genotype combination. There are 24 such cell-wise contrasts. Raw P values are retained and Benjamini–Hochberg FDR correction is applied jointly across those 24 tests. Figure significance symbols are based on the FDR-adjusted values.

### Temperature and genotype summaries

Marginal means for Temperature and Genotype are compared with Tukey adjustment and displayed as compact-letter displays using `multcomp`.

## Secondary endpoint: D54

The final disease score at 54 DAI is analysed with the same full factorial mixed model:

```text
D54 ~ Treatment * Temperature * Genotype + (1 | Rep)
```

It is treated as a secondary response, with the same cell-wise Control-vs-Treated FDR procedure.

## Repeated disease progression

For visualisation and descriptive analysis of the seven repeated disease scores, the script reshapes the data to long format and fits an `nlme::lme` model. Experimental unit is the repeated-measures subject and an AR(1) residual correlation is attempted across DAI. If the AR(1) model fails to converge, the script falls back to the same mixed model without the AR(1) term.

The script intentionally does **not** perform separate significance tests at every date, avoiding a large family of per-date tests.

## Residual diagnostics

For AUDPC and D54, the script reports Shapiro–Wilk and Anderson–Darling tests of model residuals. AUDPC residuals are additionally examined using simulated residual diagnostics from `DHARMa`. These diagnostics support evaluation of the Gaussian mixed-model assumptions but should be interpreted together with residual plots and model structure rather than as stand-alone pass/fail tests.

## Early and late disease phases

The supplied script includes a supplementary two-phase analysis intended to separate the temperature-associated early disease level from later progression:

- **Early phase:** mean disease score across 11, 18 and 26 DAI.
- **Late phase:** per-unit linear slope of disease score against DAI across 32, 40, 47 and 54 DAI.

In the supplied implementation, formal temperature tests for these two derived variables use **Control plants only**:

```text
early_level ~ Temperature * Genotype + (1 | Rep)
late_rate   ~ Temperature * Genotype + (1 | Rep)
```

A Control-only repeated-measures Temperature × Phase model is also fitted, and an AUDPC model restricted to Control plants provides the net temperature comparison.

### Manuscript-alignment issue to resolve

The current manuscript Methods describe early disease and late progression as being analysed with Treatment, Temperature, Genotype and their interactions, whereas the supplied two-phase script uses Control-only formal models. The Results text emphasizes untreated/control plants when interpreting these temperature effects. For exact reproducibility, the final manuscript and final code should state the same inferential model. This repository preserves the supplied code and flags the discrepancy rather than silently rewriting it.

## Figures and tables produced by the supplied script

Main candidate outputs include:

- genotype-specific AUDPC under each temperature, with FDR-adjusted Control-vs-Treated significance;
- global AUDPC by temperature and treatment;
- disease progression curves;
- D54 distributions;
- descriptive treatment-associated AUDPC reduction;
- a two-phase figure showing disease progression, early level, late progression rate and net AUDPC.

The script exports the full ANOVA, contrast, CLD and diagnostic tables rather than only the plotted summaries.

## Method references

- Bates D, Mächler M, Bolker B, Walker S. 2015. *Fitting Linear Mixed-Effects Models Using lme4*. Journal of Statistical Software 67:1–48. https://doi.org/10.18637/jss.v067.i01
- Kuznetsova A, Brockhoff PB, Christensen RHB. 2017. *lmerTest Package: Tests in Linear Mixed Effects Models*. Journal of Statistical Software 82:1–26. https://doi.org/10.18637/jss.v082.i13
- Benjamini Y, Hochberg Y. 1995. *Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing*. JRSS B 57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Hothorn T, Bretz F, Westfall P. 2008. *Simultaneous Inference in General Parametric Models*. Biometrical Journal 50:346–363. https://doi.org/10.1002/bimj.200810425
- Pinheiro JC, Bates DM. 2000. *Mixed-Effects Models in S and S-PLUS*. Springer. https://doi.org/10.1007/b98882
- Lenth RV et al. `emmeans`: estimated marginal means. See the repository bibliography for the software citation used by the manuscript.

The complete project bibliography is in `references/references.bib`.
