# Restricted cubic splines (3 knots) for DNAmTL, with median DNAmTL as reference.
# Default covariates reproduce the supplied spline code. That code omits
# cardiovad and diabetesT2, although the manuscript says 'same covariates'.
# To investigate the fully adjusted alternative explicitly, set
# options(dnamtl.spline_full_adjustment = TRUE) before running this script.
source("functions.R", local = TRUE)
require_packages(c("rms", "ggplot2", "ggpubr"))
library(rms)
library(ggplot2)
df_all <- complete_data(read_analysis_data(), cohort_covariates)
groups <- split_age_groups(df_all)
adjustment <- if (isTRUE(getOption("dnamtl.spline_full_adjustment", FALSE))) {
  cox_covariates
} else {
  setdiff(cox_covariates, c("cardiovad", "diabetesT2"))
}

fit_spline <- function(data, outcome) {
  # Fix knots using the original input cohort, before outcome-specific interval
  # exclusions. This preserves the source cph() knot placement for PFS.
  knots <- Hmisc::rcspline.eval(data$DNAmTL, nk = 3, knots.only = TRUE)
  spline_term <- paste0("rcs(DNAmTL, c(", paste(sprintf("%.17g", knots), collapse = ", "), "))")
  formula <- survival_formula(outcome, c(spline_term, adjustment))
  d <- model_data(data, formula, outcome)
  # Passing a datadist object avoids reliance on an object in the global workspace.
  old <- options(datadist = rms::datadist(data[, unique(c("DNAmTL", adjustment)), drop = FALSE]))
  on.exit(options(old), add = TRUE)
  fit <- rms::cph(formula, data = d, x = TRUE, y = TRUE)
  predictions <- as.data.frame(rms::Predict(fit, name = "DNAmTL", ref.zero = TRUE, fun = exp))
  tests <- as.data.frame(anova(fit))
  tests$Term <- rownames(tests)
  tests$N <- nrow(d)
  list(predictions = predictions, tests = tests)
}

plots <- list()
tests <- list()
predictions <- list()
for (outcome in c("OS", "PFS")) {
  for (group in names(groups)) {
    result <- fit_spline(groups[[group]], outcome)
    key <- paste(group, outcome, sep = "_")
    result$tests$Group <- group
    result$tests$Outcome <- outcome
    result$tests$Adjustment <- paste(adjustment, collapse = "; ")
    tests[[key]] <- result$tests
    predictions[[key]] <- transform(result$predictions, Group = group, Outcome = outcome)
    plots[[key]] <- ggplot(result$predictions, aes(DNAmTL, yhat)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.3) +
      geom_line(colour = "steelblue", linewidth = 1) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
      coord_cartesian(ylim = c(0, 4)) +
      labs(x = "DNAmTL", y = "Adjusted hazard ratio", title = paste(group, outcome)) +
      theme_classic(base_size = 12)
  }
}
suffix <- if (isTRUE(getOption("dnamtl.spline_full_adjustment", FALSE))) "full_adjustment" else "source_adjustment"
write_results(list(tests = do.call(rbind, tests), curves = do.call(rbind, predictions)),
              paste0("dose_response_", suffix, ".xlsx"))
figure <- ggpubr::ggarrange(plotlist = plots, ncol = 2, nrow = 2, labels = LETTERS[1:4])
ggsave(result_file(paste0("dose_response_", suffix, ".pdf")), figure, width = 10, height = 8)
save_session_info("11_splines")
