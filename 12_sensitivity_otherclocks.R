# Associations of alternative methylation clocks with CRC outcomes.
# Raw clock values and the supplied age-difference variables are separate analyses.
source("functions.R", local = TRUE)
df_all <- read_analysis_data()
clock_sets <- list(
  raw_clocks = c("PCHorvath1", "PCHorvath2", "PCHannum", "PCPhenoAge", "PCGrimAge",
                 "DNAmGrimAge2BasedOnRealAge", "DunedinPACE"),
  age_differences = c("diffagehorvath", "diffagehorvath2", "diffageHannum",
                      "diffagePhenoAge", "diffagepcgrim", "diffagegrim2", "DunedinPACE")
)
results <- list()
for (cutoff in c(55, 50)) {
  for (clock_set in names(clock_sets)) {
    variables <- clock_sets[[clock_set]]
    groups <- prepare_cox_groups(df_all, cutoff, variables)
    scores <- as.vector(rbind(variables, paste0(variables, "_dich")))
    result <- run_cox_grid(groups, scores)
    result$Age_cutoff <- cutoff
    result$Clock_set <- clock_set
    results[[paste(cutoff, clock_set)]] <- result
  }
}
write_results(do.call(rbind, results), "other_methylation_clocks.xlsx")
save_session_info("12_other_clocks")
