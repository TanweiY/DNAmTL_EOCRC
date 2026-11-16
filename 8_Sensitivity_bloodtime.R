# Additional adjustment for blood collection relative to surgery/chemotherapy.
source("functions.R", local = TRUE)
groups <- prepare_cox_groups(read_analysis_data())
adjustments <- c(surgery = "brelsurgtim_month", chemotherapy = "anychemorebldt")
results <- lapply(names(adjustments), function(adjustment) {
  result <- run_cox_grid(groups, outcomes = c("OS", "PFS"),
                          extra_covariates = adjustments[[adjustment]])
  result$Additional_adjustment <- adjustment
  result
})
# Chemotherapy uses a categorical timing variable, not a continuous interval.
write_results(do.call(rbind, results), "blood_collection_adjustment.xlsx")
save_session_info("08_bloodtime")
