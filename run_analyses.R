# Run the clinical analyses from a fresh session at the project root.
# Raw IDAT preprocessing is optional and intentionally run separately.
source("functions.R", local = TRUE)
# Resolve/select the dataset once for the entire run.
if (is.null(getOption("dnamtl.data_file"))) {
  options(dnamtl.data_file = input_file("df_followupdate.RData"))
}
scripts <- c(
  "2_descriptive_analyses.R", "3_robust_linear.R", "4_cox_presumption_test.R",
  "5_cox_models.R", "7_bootstrapping eocrc.R", "8_Sensitivity_bloodtime.R",
  "9_Sensitivity_Subgroup.R", "10_granular10year.R", "11_dose_response_curves.R",
  "12_sensitivity_otherclocks.R", "13_time_dependent_AUC.R", "14_auc_plots.R",
  "15_Harrell's C index.R", "interactionterms.R"
)
for (script in scripts) {
  message("Running ", script)
  source(script, local = new.env(parent = globalenv()), encoding = "UTF-8")
}
# This figure reproduces the supplied summary workbook. It is separate from
# the model outputs because its manually entered interaction P values need review.
if (file.exists(file.path("data", "forestdf.xlsx")) ||
    file.exists("forestdf.xlsx") || !is.null(getOption("dnamtl.forest_file"))) {
  source("6_forestplot_aHR.R", local = new.env(parent = globalenv()), encoding = "UTF-8")
} else {
  message("Forest figure skipped: supply forestdf.xlsx and run 6_forestplot_aHR.R separately.")
}
