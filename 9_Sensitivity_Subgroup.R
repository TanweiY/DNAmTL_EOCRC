# Exclusion analyses and clinical/molecular subgroups.
# Each restriction is followed by the original common covariate restriction,
# then age-group-specific standardization and median splitting.
source("functions.R", local = TRUE)
df_all <- read_analysis_data()
assert_columns(df_all, c("brelsurgtim_month", "anychemocat", "neotreat", "MSI_gent",
                         "brafmut", "krasmut", "cimphi", "Location"))
subsets <- list(
  exclude_2weeks_surgery = with(df_all, is.na(brelsurgtim_month) |
                                brelsurgtim_month < 0 | brelsurgtim_month > 14 / 30),
  exclude_2weeks_chemo = with(df_all, is.na(anychemocat) | anychemocat != "0-2 weeks after"),
  exclude_neoadjuvant_and_2weeks_chemo = with(df_all, neotreat == "No" &
                                             (is.na(anychemocat) | anychemocat != "0-2 weeks after")),
  MSS = df_all$MSI_gent == "No", BRAF_wildtype = df_all$brafmut == "No",
  KRAS_wildtype = df_all$krasmut == "No", CIMP_negative = df_all$cimphi == "No",
  TNM_I_II = df_all$TNMdich == "I/II", TNM_III_IV = df_all$TNMdich == "III/IV",
  rectum = df_all$Location == "Rectum", colon = df_all$Location != "Rectum",
  ever_smoker = df_all$smokingdich != "never",
  male = df_all$Sex == "Male", female = df_all$Sex == "Female"
)
results <- lapply(names(subsets), function(subgroup) {
  # which() excludes unknown subgroup membership instead of creating NA rows.
  data <- df_all[which(subsets[[subgroup]]), , drop = FALSE]
  result <- run_cox_grid(prepare_cox_groups(data))
  result$Subgroup <- subgroup
  result
})
write_results(do.call(rbind, results), "subgroup_sensitivity.xlsx")
save_session_info("09_subgroups")
