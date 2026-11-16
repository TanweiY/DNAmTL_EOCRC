# Local inputs

Place only authorized local inputs in this directory. Its contents are excluded from Git, except this README. No example participant records are provided.

## Clinical analysis file

`df_followupdate.RData` must contain a data frame named `df_all`, with one record per participant. The loader validates key columns, loads the file in an isolated environment and returns a copy. It derives age groups from `Age_diag`, preserving all original clinical values.

| Variables | Meaning / expected representation |
|---|---|
| `Age_diag`, `DNAmTL` | Numeric age at diagnosis in years and DNAmTL in the source score units |
| `time_blood` | Days from diagnosis to blood collection; delayed-entry time |
| `timeD`, `death_all` | OS exit time in days and event indicator, 0/1 |
| `timeD_PFS`, `PFS` | PFS exit time in days and progression/death indicator, 0/1 |
| `timeD_recurr`, `recurr_cp` | Recurrence exit time and cause code; recurrence is code 1 |
| `death_crccp` | Cause code for death; colorectal cancer death is code 1 |
| `timey`, `timey_PFS` | Corresponding OS/PFS times from diagnosis in years, used for prediction evaluation |
| `Sex` | Factor, including `Female` and `Male` |
| `TNMdich`, `TNM_adj` | Factors for I/II versus III/IV, and separate I-IV stages |
| `locationdich`, `Location` | Source proximal/distal grouping and detailed tumour location, including `Proximal colon`, `Distal colon`, `Rectum` |
| `chemradther` | Source adjuvant-treatment factor |
| `cardiovad`, `diabetesT2`, `high_blood_pressure`, `other_cancer` | Source medical-history factors |
| `Schooling_years`, `Schooling_yearsdich` | Source education categories; dichotomy is <9 versus >=9 years, not <=9 versus >9 |
| `smoking`, `smokingdich` | Detailed smoking status and `never`/`ever` grouping |
| `lifetimealcohol_day`, `lifealcoholdich`, `liftimdrink_cat` | Continuous and grouped lifetime alcohol measures |
| `BMI`, `BMI_cat`, `BMI_5_14earlier`, `BMI5_14er_cat` | Source current and prediagnostic BMI measures/categories |
| `brelsurgtim_month` | Blood-collection time relative to surgery, in months |
| `brelsurgtim_monthcat`, `anychemocat` | Timing categories; chemotherapy exclusion uses `0-2 weeks after` |
| `anychemorebldt` | Categorical chemotherapy timing: `afteranychemo`, `beforeanychemo`, `noanychemo` |
| `neotreat` | Neoadjuvant treatment, including `No` and `Yes` |
| `MSI_gent`, `brafmut`, `krasmut`, `cimphi` | Molecular subgroup factors; negative group is `No` |
| `TL_shortenmeanbirh` | Supplied descriptive DNAmTL shortening index; retained rather than recalculated |
| `tn_id` | Participant identifier used only in the original robust-regression eligibility restriction; never exported |

The alternative-clock script additionally requires `PCHorvath1`, `PCHorvath2`, `PCHannum`, `PCPhenoAge`, `PCGrimAge`, `DNAmGrimAge2BasedOnRealAge`, `DunedinPACE`, `diffagehorvath`, `diffagehorvath2`, `diffageHannum`, `diffagePhenoAge`, `diffagepcgrim` and `diffagegrim2`.

Factors and reference levels are retained from the supplied data, except `Location` is explicitly referenced to `Proximal colon` in robust regression. Do not replace unknown values with zero or silently recode factors to integers.

## Figure inputs

`forestdf.xlsx` is a summary workbook with columns `Group`, `N`, `HR`, `LL`, `UL` and `P for interaction`. Its heading rows are `Overall mortality`, `Disease progression or death`, `55-year cut-off` and `50-year cut-off`; estimate rows use `EOCRC`/`LOCRC`. Review manually entered sample sizes and interaction values against current model outputs before using the figure for publication.

`A4_portrait_blank_template.pptx` is optional. It must contain a Blank layout and have a 210 x 297 mm slide size.

## Raw methylation inputs

Script 1 expects IDAT files and `SampleSheet.csv` under the project-root `rawdata/` directory. The sample sheet must identify samples, Sentrix IDs/positions, batch `Pool_ID` and age. These files are not included. ChAMP must be installed separately. Preprocessing has substantially greater memory and runtime needs than the clinical analyses.
