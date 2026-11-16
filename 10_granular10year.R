# Event counts and adjusted associations across finer age bands.
source("functions.R", local = TRUE)
require_packages("dplyr")
df_all <- read_analysis_data()
df_all$age_band <- cut(df_all$Age_diag, c(-Inf, 30, 40, 50, 60, 70, 80, Inf),
                       labels = c("<30", "30-39", "40-49", "50-59", "60-69", "70-79", ">=80"),
                       right = FALSE)
counts <- dplyr::summarise(dplyr::group_by(df_all, age_band),
                          N = dplyr::n(), OS_events = sum(death_all == 1, na.rm = TRUE),
                          PFS_events = sum(PFS == 1, na.rm = TRUE),
                          OS_missing = sum(is.na(death_all)), PFS_missing = sum(is.na(PFS)), .groups = "drop")
write_results(counts, "age10_event_summary.xlsx")
df_all$age_band2 <- cut(df_all$Age_diag, c(-Inf, 50, 60, 70, 80, Inf),
                        labels = c("<50", "50-59", "60-69", "70-79", ">=80"), right = FALSE)
data <- complete_data(df_all, cohort_covariates)
groups <- lapply(split(data, data$age_band2, drop = TRUE), standardize_scores)
write_results(run_cox_grid(groups, scores = "DNAmTL", outcomes = c("OS", "PFS")),
              "aHR_age10band.xlsx")

# Report distinct questions separately: categorical age-band interaction,
# with/without within-band age adjustment, and continuous age interaction.
tests <- list()
for (outcome in c("OS", "PFS")) {
  for (adjust_age in c(TRUE, FALSE)) {
    vars <- if (adjust_age) cox_covariates else setdiff(cox_covariates, "Age_diag")
    result <- interaction_lrt(data, "age_band2", outcome, covariates = vars)
    result$Adjustment <- if (adjust_age) "Includes continuous age" else "Age bands only"
    tests[[paste(outcome, adjust_age)]] <- result
  }
  result <- interaction_lrt(data, "Age_diag", outcome)
  result$Adjustment <- "Continuous age interaction"
  tests[[paste(outcome, "continuous")]] <- result
}
write_results(do.call(rbind, tests), "age_interaction_tests.xlsx")
save_session_info("10_age_bands")
