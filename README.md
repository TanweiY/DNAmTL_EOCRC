# DNAmTL_EOCRC

R code for studying DNA methylation-based telomere length (DNAmTL) and colorectal cancer prognosis in the DACHS cohort. The primary age-of-onset comparison uses <55 versus >=55 years; sensitivity analyses use <50 versus >=50 years.

This repository contains analysis code and documentation. Participant-level data, raw methylation arrays, manuscripts, and generated results are not distributed.

## Quick start

1. Open `DNAmTL.Rproj` in RStudio, or start R in this repository's root directory.
2. Run `source("install_dependencies.R")` once if required packages are missing.
3. Put the authorized analysis file `df_followupdate.RData` in `data/`. It must contain a data frame named `df_all`. In an interactive R session, a file picker is offered when the file is absent; no personal directory needs to be written into the code.
4. Run one analysis, for example:

```r
source("5_cox_models.R", encoding = "UTF-8")
```

Each analysis sources the shared `functions.R` and loads a fresh copy of the data. Outputs are written to `results/`, with session information identifying the R and package versions used. Input RData files are never overwritten.

To run all clinical analyses:

```r
source("run_analyses.R", encoding = "UTF-8")
```

The forest plot additionally needs the supplied summary workbook `data/forestdf.xlsx`. An optional `data/A4_portrait_blank_template.pptx` gives the AUC PowerPoint a true A4 portrait slide; otherwise the figure is fitted to a standard slide. The PDF remains A4.

Batch execution from the project root also works:

```sh
Rscript validate_code.R
Rscript run_analyses.R
```

For authorized data stored outside the project, an in-session option can override the default. This is optional, never hard-coded into repository files:

```r
options(dnamtl.data_file = file.choose())
source("run_analyses.R", encoding = "UTF-8")
```

## Script guide

| File | Analysis |
|---|---|
| `1_DNAmethylation_preprocessing.R` | Optional EPIC-array import, filtering, BMIQ normalization and ComBat |
| `2_descriptive_analyses.R` | Baseline/timing tables, correlations, scatterplots, shortening-index comparisons and slope interaction |
| `3_robust_linear.R` | Robust regression and Type II global Wald tests |
| `4_cox_presumption_test.R` | Schoenfeld diagnostics, log-minus-log plots and time-varying DNAmTL effects |
| `5_cox_models.R` | Delayed-entry Cox models at both age cutoffs, reverse-KM follow-up, descriptive mortality curves |
| `6_forestplot_aHR.R` | Forest plot from an explicitly supplied summary workbook |
| `7_bootstrapping eocrc.R` | 1,000 EOCRC bootstrap resamples for each of OS and PFS |
| `8_Sensitivity_bloodtime.R` | Additional adjustment for surgery/chemotherapy timing |
| `9_Sensitivity_Subgroup.R` | Timing exclusions and clinical/molecular subgroups |
| `10_granular10year.R` | Age-band counts, Cox estimates and distinct age interaction tests |
| `11_dose_response_curves.R` | Three-knot restricted cubic spline curves and nonlinearity tests |
| `12_sensitivity_otherclocks.R` | Other methylation clocks and supplied age-difference measures |
| `13_time_dependent_AUC.R` | Paired apparent AUCs at both age cutoffs |
| `14_auc_plots.R` | Eight-panel AUC figure; run script 13 first |
| `15_Harrell's C index.R` | Paired apparent C-indices and comparisons |
| `interactionterms.R` | Age-group interactions, explicitly separating pooled and within-group scaling |
| `functions.R` | Shared loading, validation, transformations, model fitting and result extraction |

## Analysis conventions

- OS and PFS are the main outcomes. CSS and TTR are retained as supplementary analyses from the supplied code. CSS/TTR models use cause-specific event indicators (`== 1`); they are not Fine-Gray models.
- Association models use time since diagnosis in days, with entry at blood collection (`time_blood`). Invalid entry/exit intervals and missing model variables are excluded and reported through fitted sample sizes.
- The common Cox covariate restriction is applied before age-group-specific standardization, matching the supplied analyses. Continuous exposures represent one within-group SD. Median categories are exploratory.
- `N_input`, fitted `N`, and `N_events` are separate. `Fit_status`, `Fit_warnings`, and `omitted_covariates` identify numerical warnings and covariates constant within a subgroup. A completed script does not establish that every fitted model is statistically reliable.
- Robust-regression global P values are Type II Wald tests. Numeric estimates, confidence limits and P values are exported alongside display strings. Model sizes are calculated rather than entered manually.
- AUC and C-index analyses retain the supplied right-censored, same-data evaluation design. They do **not** implement delayed entry, external validation or optimism correction. Bootstrap coefficient stability is a different analysis from prediction-model validation.
- The source spline model omits cardiovascular disease and diabetes, unlike the main Cox model. This specification is retained explicitly for reproducibility. To examine the alternative with all primary covariates, run script 11 after `options(dnamtl.spline_full_adjustment = TRUE)`; its outputs use a different suffix.
- Interaction P values depend on the exposure scale and interaction definition. Continuous, median-categorized, categorical age-band and continuous-age tests are exported separately. The forest workbook is not automatically used as statistical evidence.

See [data/README.md](data/README.md) for input definitions and [REVIEW_NOTES.md](REVIEW_NOTES.md) for changes and unresolved methodological questions.

## Scope and verification

The input RData already contains derived clinical variables, DNAmTL, other clocks and outcomes. The original files did not include the full pipeline that creates this analysis dataset or computes clocks from methylation. Script 1 therefore cannot, on its own, generate the input needed by scripts 2-15.

`validate_code.R` checks syntax and key behavior with synthetic data. Full clinical validation additionally requires authorized local data. See [VALIDATION.md](VALIDATION.md) for the validation performed on this revision.

Package documentation: [survival](https://CRAN.R-project.org/package=survival), [riskRegression](https://CRAN.R-project.org/package=riskRegression), [compareC](https://CRAN.R-project.org/package=compareC), [rms](https://CRAN.R-project.org/package=rms), [ChAMP](https://bioconductor.org/packages/ChAMP/).
