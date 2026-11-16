# Paired apparent Harrell C-indices for the same two prediction models as AUC.
# Right-censored times from diagnosis are retained from the supplied analysis.
source("functions.R", local = TRUE)
df_all <- read_analysis_data()
outcomes <- list(OS = c(time = "timeD", event = "death_all"),
                  PFS = c(time = "timeD_PFS", event = "PFS"))
results <- list()
for (cutoff in c(55, 50)) {
  groups <- split_age_groups(df_all, cutoff)
  for (group in names(groups)) {
    for (outcome in names(outcomes)) {
      spec <- outcomes[[outcome]]
      label <- paste(group, if (group == "EOCRC") "<" else ">=", cutoff)
      result <- get_cindex_results(groups[[group]], spec[["time"]], spec[["event"]], label, outcome)
      result$Age_cutoff <- cutoff
      results[[paste(cutoff, group, outcome)]] <- result
    }
  }
}
write_results(do.call(rbind, results), "Cindex_comp.xlsx")
save_session_info("15_Cindex")
