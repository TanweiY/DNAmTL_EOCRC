# Associations of baseline characteristics with standardized DNAmTL.
# Each age group has its own SD; global factor P values are Type II Wald tests.
source("functions.R", local = TRUE)
df_all <- read_analysis_data()
df_all$Location <- relevel(factor(df_all$Location), ref = "Proximal colon")

x_vars <- c("Age_diag", "Sex", "Schooling_years", "BMI_5_14earlier", "smoking",
            "lifealcoholdich", "cardiovad", "diabetesT2", "high_blood_pressure",
            "other_cancer", "TNM_adj", "Location")
# Retain the original pre-standardization restriction. The helper then excludes
# missing outcomes/other_cancer for each fitted model and reports the actual N.
scaling_vars <- c("tn_id", setdiff(x_vars, "other_cancer"))
df_complete <- complete_data(df_all, scaling_vars)
groups <- lapply(split_age_groups(df_complete), standardize_scores)
results <- lapply(names(groups), function(group) {
  result <- multivariable_rlm_single(groups[[group]], c(DNAmTL = "DNAmTL"), x_vars)
  result$Group <- group
  result$Global_P_display <- format_p(result$global_p_value)
  result
})
long <- do.call(rbind, results)
# Join by coefficient names; factors need not have identical rows in both groups.
display <- lapply(split(long, long$Group), function(d) {
  d[, c("Y_Variable", "Term", "N", "estimateCI", "Global_P_display")]
})
wide <- merge(display$EOCRC, display$LOCRC, by = c("Y_Variable", "Term"),
               all = TRUE, suffixes = c("_EOCRC", "_LOCRC"), sort = FALSE)
write_results(list(display = wide, numeric_results = long), "robustl_globalp.xlsx")
save_session_info("03_robust_linear")
