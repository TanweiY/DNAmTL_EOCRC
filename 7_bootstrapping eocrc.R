# EOCRC coefficient stability: 1,000 participant-level bootstrap resamples.
# This describes sampling variability; it is not out-of-sample validation or an
# optimism-corrected estimate of predictive performance.
source("functions.R", local = TRUE)
require_packages(c("survival", "ggplot2", "ggpubr"))
library(ggplot2)
eo <- standardize_scores(split_age_groups(read_analysis_data())$EOCRC)
B <- getOption("dnamtl.bootstrap_replicates", 1000L)
stopifnot(length(B) == 1, is.finite(B), B >= 2, B == as.integer(B))

boot_results <- list()
summaries <- list()
plots <- list()
for (outcome in c("OS", "PFS")) {
  formula <- survival_formula(outcome, c("DNAmTL", cox_covariates))
  # Match the original resampling population: resample the full EOCRC group;
  # coxph performs complete-case exclusion within each bootstrap sample.
  set.seed(123)
  samples <- lapply(seq_len(B), function(i) {
    d <- eo[sample.int(nrow(eo), nrow(eo), replace = TRUE), , drop = FALSE]
    warnings <- character()
    fit <- tryCatch(withCallingHandlers(
      survival::coxph(formula, data = d),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }
    ), error = function(e) e)
    error <- if (inherits(fit, "error")) conditionMessage(fit) else ""
    beta <- if (inherits(fit, "error")) NA_real_ else unname(coef(fit)["DNAmTL"])
    # Retain every replicate and its warning/error status for inspection.
    data.frame(Outcome = outcome, Replicate = i, Log_HR = beta, HR = exp(beta),
               Warning = paste(warnings, collapse = "; "), Error = error)
  })
  result <- do.call(rbind, samples)
  # Preserve the original summary over finite estimates, and expose warnings.
  # Filtering out every warning would discard all PFS replicates merely because
  # Surv() reports the known invalid intervals that it then excludes internally.
  valid <- is.finite(result$HR) & result$Error == ""
  numerical_warning <- grepl("converg|infinite|iteration|singular", result$Warning, ignore.case = TRUE)
  result$Numerical_warning <- numerical_warning
  if (!any(valid)) stop("No successful bootstrap fits for ", outcome)
  hr <- result$HR[valid]
  q <- quantile(hr, c(0.025, 0.25, 0.5, 0.75, 0.975))
  summaries[[outcome]] <- data.frame(
    Outcome = outcome, B_requested = B, B_finite = sum(valid),
    B_failed_or_nonfinite = sum(!valid), B_numerical_warning = sum(numerical_warning),
    Mean_HR = mean(hr), Median_HR = q[3],
    Q25 = q[2], Q75 = q[4], Percentile_CI_lower = q[1], Percentile_CI_upper = q[5]
  )
  boot_results[[outcome]] <- result
  plots[[outcome]] <- ggplot(data.frame(HR = hr), aes(HR)) +
    geom_histogram(bins = 30, colour = "black", fill = "lightblue") +
    geom_vline(xintercept = q[2:4], linetype = "dashed", colour = "red") +
    labs(title = if (outcome == "OS") "Overall mortality" else "Progression or death",
         subtitle = sprintf("Median %.2f; IQR %.2f-%.2f; %d/%d finite fits (%d numerical warnings)",
                            q[3], q[2], q[4], sum(valid), B, sum(numerical_warning)),
         x = "Adjusted hazard ratio per SD", y = "Count") + theme_classic()
}
write_results(list(summary = do.call(rbind, summaries), replicates = do.call(rbind, boot_results)),
              "bootstrap_EOCRC.xlsx")
figure <- ggpubr::ggarrange(plotlist = plots, ncol = 1, labels = c("A", "B"))
ggsave(result_file("bootstrap_EOCRC.pdf"), figure, width = 8, height = 10)
save_session_info("07_bootstrap")
