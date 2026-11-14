# Transfer Learning for High-Dimensional Generalized Linear Models

This project studies Trans-GLM using a Gaussian simulation and prepares a gene-expression dataset for a future gemcitabine-response analysis. The simulation uses six candidate sources, 0–6 transferable sources, `h = 1, 10, 20`, and 10 repetitions per setting. The real-data preprocessing has now been completed; fitting Trans-GLM on the real data is a separate next step.

## Repository files

| File | Purpose |
|---|---|
| `00_requirements.R` | Loads the R packages needed for simulation, data processing, and reports. |
| `01_SimulationFunctions.R` | Functions used to generate simulated observations and calculate MSE and MAE. |
| `02_GaussianSimulation.R` | Fits the three simulation methods and saves errors and figures. |
| `03_PrepareGDSC2.R` | Rebuilds the cleaned gene-expression dataset from the three original Cell Model Passports downloads. |
| `10_Simulation_DataProcessing.Rmd` | Demonstrates simulated data generation and plots the saved simulation results. |
| `20_RealData_DataProcessing.Rmd` | Documents the real-data merge and gene selection and checks the cleaned data. |
| `FinalReport.Rmd` | R Markdown report; the font dependency was removed for easier PDF knitting. |
| `FinalReport.pdf` | PDF report. |
| `data/gdsc2_gemcitabine_2000genes.csv` | Final modeling table: 229 cell models, `LN_IC50`, and 2,000 log2(TPM+1) predictors. |
| `data/gdsc2_gemcitabine_samples.csv` | Matched model IDs, group labels, cancer types, and outcomes. |
| `data/gdsc2_gemcitabine_merge_audit.csv` | All 233 eligible drug-response models and whether a matching TPM column exists. |
| `data/gdsc2_merge_counts.csv` | Number of drug-response and matched-expression models by group. |
| `data/gdsc2_gene_filter_counts.csv` | Numbers of genes retained during filtering. |
| `data/gdsc2_selected_genes.csv` | Names, identifiers, and source-based variance rankings of the selected 2,000 genes. |
| `results/` | The four CSV output files from the preliminary Gaussian simulation. |
| `figures/` | Two figures from the preliminary Gaussian simulation. |

## Packages

Install once in R if necessary:

```r
install.packages(c("glmtrans", "glmnet", "ggplot2", "rmarkdown", "readxl", "data.table"))
```

Then set the working directory to the repository folder in RStudio and run:

```r
source("00_requirements.R")
```

The `data.table` package is used to read just the matching columns of the large TPM file, rather than loading all cell models.

## Reproduce the Gaussian simulation

1. In RStudio, set the working directory to the top-level repository folder.
2. Run `source("02_GaussianSimulation.R")` to regenerate results and figures. The script loads the packages and helper functions itself.
3. Knit `10_Simulation_DataProcessing.Rmd` to check the resulting summaries and plots.
4. Knit `FinalReport.Rmd` to PDF if a LaTeX distribution (e.g. TinyTeX) is installed.

The simulation uses 10 repetitions for each of 21 combinations of `h` and `Ka`. There are 6 source datasets of 100 observations each, 100 target training observations, 1,000 target test observations, 500 predictors, and 10 nonzero coefficients. Three models are compared: target-only Lasso, Trans-GLM with all sources, and Trans-GLM with automatic source screening. The seed formulas depend on repetition and source count, permitting reproducible draws and a controlled comparison across values of `h`.

## Reproduce the GDSC2 gene-expression preprocessing

The repository contains the cleaned dataset but **not the original large downloads**. To independently regenerate the cleaned CSVs, download the following public inputs from [Cell Model Passports](https://cellmodelpassports.sanger.ac.uk/downloads/gdsc) and put them in a local folder named `data/raw/`:

- `GDSC2_fitted_dose_response_27Oct23.xlsx`
- `model_list_20261013.csv`
- `rnaseq_merged_rsem_tpm_20260323.csv` (the TPM version, **not** FPKM or expected counts)

Then run:

```r
source("03_PrepareGDSC2.R")
```

The script filters GDSC2 to gemcitabine, matches `SANGER_MODEL_ID`/`model_id`, retains pancreatic, colorectal, gastric, non-small-cell lung and ovarian cancers, and drops models without RNA-seq TPM. It then removes genes with invalid or missing expression, retains genes with TPM > 1 in at least 20% of **external-source** samples, and takes the 2,000 highest-variance genes based on **source-only** log2(TPM+1) expression. The resulting data have 229 unique models and 2,000 predictors. The target sample was not used to select genes. These predictors are **not yet standardized** for model fitting; centering/scaling and any target-based tuning must occur inside training folds.

Knit `20_RealData_DataProcessing.Rmd` to examine the saved counts and cleaned data. This notebook can be knitted from the repository without downloading the original files again. The raw-input script can be run separately when the original files are available.

### Real-data sample counts

| Cancer group | Trans-GLM role | Drug-response rows | Matched with TPM |
|---|---|---:|---:|
| Pancreatic | Target | 30 | 30 |
| Colorectal | Source | 47 | 47 |
| Gastric | Source | 28 | 28 |
| Non-small-cell lung | Source | 88 | 84 |
| Ovarian | Source | 40 | 40 |
| **Total** | | **233** | **229** |

Four lung models lack matching TPM measurements and were excluded. The TPM file originally contains 41,145 gene rows; 36,417 have valid values in all matched models, 15,933 pass the source-expression threshold, and the highest-variance 2,000 are retained.

## GitHub upload guidance

**Upload the repository contents**, including the two R Markdown notebooks, all `.R` scripts, the report files, `figures/`, `results/`, and the six small prepared data/metadata CSVs. The modeling dataset is about 6 MB, so it fits normal GitHub file limits.

**Do not upload the raw 355 MB RNA-seq TPM file, the FPKM/count alternatives, or the original 22 MB GDSC2 workbook** to this repository. These original downloads remain on your computer under `data/raw/`, which is ignored by `.gitignore`. The processing script and the input-file names document how to obtain and rebuild the smaller data.

## Current status and limitations

The preliminary Gaussian simulation and initial GDSC2 cleaning are complete, but the gene-expression data have not yet been fitted using Trans-GLM. No real-data prediction accuracy, selected source groups, or biological findings are claimed here. Before reporting final performance, use held-out pancreatic models and keep target-data preprocessing and hyperparameter selection inside the training splits. Independently verify the script output with a fresh local R run; R execution could not be tested in the environment where the cleaned CSVs were prepared.

## References

Tian, Y., & Feng, Y. (2023). *Transfer learning under high-dimensional generalized linear models*. Journal of the American Statistical Association. https://doi.org/10.1080/01621459.2022.2071278

Craig, E., et al. (2025). *Pretraining and the lasso*. Journal of the Royal Statistical Society Series B. https://doi.org/10.1093/jrsssb/qkaf050
