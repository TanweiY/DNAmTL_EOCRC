# Run once, explicitly, before the analyses. Existing packages are retained.
# DNA methylation preprocessing dependencies are optional and installed separately.
packages <- c(
  "survival", "MASS", "car", "tableone", "dplyr", "tidyr", "readxl", "writexl",
  "ggplot2", "ggpubr", "ggpmisc", "forestplot", "officer", "rvg",
  "rms", "Hmisc", "riskRegression", "compareC", "prodlim"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")

# Only if raw EPIC-array preprocessing is required:
# install.packages("BiocManager")
# BiocManager::install("ChAMP")
