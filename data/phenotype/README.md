# Phenotype data

This directory contains the phenotype dataset supplied with the analysis.

## Files

- `Data_phenotype.xlsx` — original workbook used by `scripts/phenotype/01_phenotype_analysis.R`.
- `Data_phenotype.csv` — Git-readable mirror of the same worksheet for inspection and version control.

## Experimental structure

The disease dataset is a complete balanced factorial design with **192 experimental units**:

- 2 treatments: `Water` and `PRI` (recoded in the script as `Control` and `Treated`);
- 4 thermal regimes: `T0`, `T1`, `T2`, `T3` (recoded as `NH`, `HS-7`, `HS-2`, `HS`);
- 6 genotypes: `DEEP`, `NANT`, `NEVA`, `OXHELLA`, `PRESTO`, `ROBILA`;
- 4 block/replicate levels;
- 7 disease assessments: 11, 18, 26, 32, 40, 47 and 54 DAI.

There are four experimental units in every Treatment × Temperature × Genotype cell and no missing disease-score values in the supplied workbook.

The workbook also contains `Root weight`, `Fresh weight` and `Dry weight`. The supplied disease-analysis script reads these columns during import but then drops them before the disease modelling shown in the manuscript; they are therefore **not part of the phenotype inference documented here**.
