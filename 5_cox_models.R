# Delayed-entry Cox models: DNAmTL per within-group SD and above versus at/below
# the within-group median. Primary cutoff 55 years; sensitivity cutoff 50 years.
source("functions.R", local = TRUE)
require_packages(c("survival", "prodlim"))
df_all <- read_analysis_data()

results <- lapply(c(55, 50), function(cutoff) {
  groups <- prepare_cox_groups(df_all, cutoff)
  result <- run_cox_grid(groups)
  result$Age_cutoff <- cutoff
  result
})
cox_results <- do.call(rbind, results)
write_results(cox_results, "multicox_all_results.xlsx")
write.csv(cox_results, result_file("multicox_all_results.csv"), row.names = FALSE)

# Reverse Kaplan-Meier follow-up estimates reproduce the supplied descriptive
# analysis from diagnosis. These are distinct from delayed-entry association models.
groups <- c(list(All = df_all), split_age_groups(df_all))
followup <- lapply(groups, function(d) {
  prodlim::prodlim(prodlim::Hist(timey, death_all) ~ 1, data = d, reverse = TRUE)
})
capture.output(lapply(followup, quantile), file = result_file("reverse_KM_followup.txt"))

# Unadjusted descriptive cumulative mortality curves from diagnosis.
# These reproduce the original visualization and do not implement delayed entry.
groups <- lapply(split_age_groups(df_all), dichotomize_scores)
for (group in names(groups)) groups[[group]]$plot_group <- paste(group, groups[[group]]$DNAmTL_dich)
plot_data <- do.call(rbind, groups)
km <- survival::survfit(survival::Surv(timey, death_all) ~ plot_group, data = plot_data)
pdf(result_file("KM_four_groups.pdf"), width = 9, height = 6)
plot(km, fun = "event", col = c("#E41A1C", "#A50F15", "#377EB8", "#08519C"),
     lty = c(1, 2, 1, 2), xlim = c(0, 15), xlab = "Time since diagnosis (years)",
     ylab = "Cumulative all-cause mortality probability")
legend("topleft", sub("plot_group=", "", names(km$strata)),
       col = c("#E41A1C", "#A50F15", "#377EB8", "#08519C"), lty = c(1, 2, 1, 2), bty = "n")
dev.off()
save_session_info("05_cox")
