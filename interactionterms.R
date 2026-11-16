# Likelihood-ratio tests of DNAmTL-by-age-group interaction.
# Pooled-SD and within-group-SD interactions test different contrasts; report
# them separately. This prevents accidental reuse of 55-year scaling for a
# 50-year interaction, which occurred in the supplied exploratory script.
source("functions.R", local = TRUE)
df_all <- read_analysis_data()
results <- list()
for (cutoff in c(55, 50)) {
  for (scaling in c("pooled", "within_age_group")) {
    if (scaling == "pooled") {
      data <- standardize_scores(dichotomize_scores(df_all))
    } else {
      data <- do.call(rbind, lapply(split_age_groups(df_all, cutoff), function(d) {
        standardize_scores(dichotomize_scores(d))
      }))
    }
    data$age_group <- age_group(data$Age_diag, cutoff)
    for (outcome in c("OS", "PFS")) {
      for (score in c("DNAmTL", "DNAmTL_dich")) {
        result <- interaction_lrt(data, "age_group", outcome, score)
        result$Age_cutoff <- cutoff
        result$Scaling <- scaling
        results[[paste(cutoff, scaling, outcome, score)]] <- result
      }
    }
  }
}
write_results(do.call(rbind, results), "DNAmTL_age_interactions.xlsx")
save_session_info("interaction")
