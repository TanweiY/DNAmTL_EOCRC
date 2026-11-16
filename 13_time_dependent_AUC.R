# Apparent time-dependent AUC: age + sex + TNM stage, with/without DNAmTL.
# Both models use the same complete cases. These right-censored prediction
# analyses reproduce the supplied design and do not account for delayed entry.
# They are not external validation or optimism-corrected performance estimates.
source("functions.R", local = TRUE)
require_packages(c("survival", "riskRegression", "dplyr", "tidyr"))
library(survival)
df_all <- read_analysis_data()
outcomes <- list(OS = c(time = "timey", event = "death_all"),
                  PFS = c(time = "timey_PFS", event = "PFS"))
auc_results <- list()
contrast_results <- list()
for (cutoff in c(55, 50)) {
  groups <- split_age_groups(df_all, cutoff)
  for (group in names(groups)) {
    for (outcome in names(outcomes)) {
      spec <- outcomes[[outcome]]
      d <- complete_data(groups[[group]], c(spec, "Age_diag", "Sex", "TNM_adj", "DNAmTL"))
      # Score's formula parser requires the unqualified Surv constructor.
      response <- paste0("Surv(", spec[["time"]], ", ", spec[["event"]], ")")
      base_formula <- as.formula(paste(response, "~ TNM_adj + Age_diag + Sex"))
      extended_formula <- as.formula(paste(response, "~ TNM_adj + DNAmTL + Age_diag + Sex"))
      models <- list(ref = survival::coxph(base_formula, data = d, x = TRUE, y = TRUE),
                     mtl = survival::coxph(extended_formula, data = d, x = TRUE, y = TRUE))
      score <- riskRegression::Score(models, formula = as.formula(paste(response, "~ 1")),
                                      data = d, conf.int = TRUE, metrics = "AUC",
                                      times = seq(0.5, 20, by = 0.5), contrasts = TRUE,
                                      split.method = "none")
      auc <- as.data.frame(score$AUC$score)
      contrast <- as.data.frame(score$AUC$contrasts)
      for (field in c("group", "outcome", "Age_cutoff", "N")) {
        value <- switch(field, group = group, outcome = outcome, Age_cutoff = cutoff, N = nrow(d))
        auc[[field]] <- value
        contrast[[field]] <- value
      }
      key <- paste(cutoff, group, outcome, sep = "_")
      auc_results[[key]] <- auc
      contrast_results[[key]] <- contrast
      directory <- paste0("auc_cutoff", cutoff)
      write_results(auc, file.path(directory, paste0("AUC", tolower(group), "_", tolower(outcome), ".xlsx")))
      write_results(contrast, file.path(directory, paste0("diffAUC", tolower(group), "_", tolower(outcome), ".xlsx")))
    }
  }
}
all_auc <- do.call(rbind, auc_results)
all_contrasts <- do.call(rbind, contrast_results)
write_results(list(AUC = all_auc, contrasts = all_contrasts), "AUC_all_results.xlsx")
annual <- subset(all_auc, times >= 1 & times <= 10 & times %% 1 == 0)
annual$AUC_CI <- sprintf("%.2f (%.2f, %.2f)", annual$AUC, annual$lower, annual$upper)
wide <- tidyr::pivot_wider(annual[, c("Age_cutoff", "outcome", "times", "group", "model", "AUC_CI")],
                            names_from = c("model", "group"), values_from = "AUC_CI")
write_results(wide, "AUC_allsuptable.xlsx")
save_session_info("13_AUC")
