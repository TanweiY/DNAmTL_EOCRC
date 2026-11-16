# Code revision notes

## Readability and execution

- Replaced Chinese and exploratory comments with English explanations of inputs, methods, output and interpretation.
- Removed personal Windows/server paths, analysis-time package installation and writes back to the clinical input file.
- Every analysis sources the same local `functions.R`. The missing chemotherapy-adjustment helper was recovered from the supplied external Cox-function source.
- Consolidated repeated model loops and transformations; retained the original filenames for traceability.
- Removed unused exploratory optimal-cutpoint bootstrap helpers. They were not called by the supplied analysis scripts and do not form part of the updated manuscript's coefficient-bootstrap analysis.
- Kept full-precision numeric outputs; display formatting is separate. Fixed P-value formatting that previously rounded first and could incorrectly label a value as `<0.0001`.
- Replaced row-position joins with term-key joins and reported actual fitted sample sizes.
- Retained warning text in model and bootstrap exports. Constant subgroup-defining covariates are omitted explicitly and listed in the output.

## Corrections with potential analytical consequences

- Corrected the alternative-clock 50-year comparison from `<50` versus `>=55` to `<50` versus `>=50`. The corrected older group includes ages 50-54.
- Made both AUC and C-index models use the same outcome-specific complete cases. AUC now consistently uses `Age_diag`; `Age` and `Age_diag` were identical in the locally validated dataset.
- Negated Cox linear predictors for `compareC`, so its concordance direction agrees with `survival::concordance`. The paired comparison P values are unchanged by negating both scores.
- Corrected the descriptive subgroup boxplot's test grouping to the displayed two age groups and the Kruskal-Wallis test described in the manuscript. The original helper used a different three-group variable and Wilcoxon test. Plot limits no longer remove observations before computing tests.
- Fixed duplicate sensitivity output names so separate exclusions no longer overwrite one another.
- Raw preprocessing aligns metadata by retained sample IDs instead of deleting a hard-coded sample. Age 50 is assigned to the intermediate age group. These upstream changes require verification with raw data before regenerating methylation inputs.

## Items requiring scientific reconciliation

The code does not silently rewrite analysis choices to match a reported number. The local author review report identifies numerical examples requiring reconciliation:

1. **Fitted sample sizes and event counts.** The source scripts used input row counts in some tables. Invalid PFS/recurrence delayed-entry intervals and additional missing predictors reduce the fitted populations. The revised exports distinguish these quantities.
2. **Interaction definitions.** Continuous and dichotomized DNAmTL interaction results, and categorical versus continuous age interactions, were mixed in exploratory code and reported materials. Use the output corresponding to the intended contrast.
3. **Adjustment sets.** Primary code includes cardiovascular disease and diabetes. The spline source omits them. The manuscript's descriptions should reflect the model actually used. A separately labeled fully adjusted spline option is available.
4. **Education threshold.** The source data/code uses `<9` versus `>=9` years. This differs from `<=9` versus `>9` wording.
5. **Treatment timing.** Chemotherapy adjustment is categorical, while surgery timing is continuous in months. Avoid describing both as continuous time intervals.
6. **Prediction evaluation.** AUC and C-index are apparent performance from right-censored models fit and evaluated on the same data; they do not use delayed entry. This differs from the primary association models. No external-validation claim is supported by these scripts.
7. **Numerical warnings.** Some small-subgroup or supplementary Cox fits generate convergence/separation warnings. Their exported coefficients should be reviewed before interpretation; this refactor does not choose penalization or alter the scientific adjustment set automatically.
8. **Incomplete upstream pipeline.** Clock calculation and assembly of the clinical analysis file are not supplied. Raw-array preprocessing is documented but was not executed in this revision.

Preserved conventions include the original common covariate restriction before SD calculation, seed 123 for each outcome's bootstrap, 1,000 replicates, source spline adjustment, source-population spline knots, and right-censored prediction evaluation. These choices are explicit so any later methodological revision can be distinguished from code cleanup.
