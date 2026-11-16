###############################################################################
# DNAmTL: Proportional hazards diagnostics (A4, eight panels)
# Panels A/C/E/G: Schoenfeld residual plots
# Panels B/D/F/H: Log-minus-log survival plots
###############################################################################

library(survival)

source("functions.R", local = TRUE)
df_all <- read_analysis_data()

adjust_vars <- c(
  "Age_diag", "Sex", "Schooling_yearsdich", "TNMdich",
  "locationdich", "chemradther", "cardiovad", "diabetesT2",
  "smokingdich", "lifetimealcohol_day"
)

eocrc <- subset(df_all, age_group55 == "EOCRC")
locrc <- subset(df_all, age_group55 == "LOCRC")

# Retain original subgroup-specific standardisation
eocrc$DNAmTL <- as.numeric(scale(eocrc$DNAmTL))
locrc$DNAmTL <- as.numeric(scale(locrc$DNAmTL))

data_groups <- list(EOCRC = eocrc, LOCRC = locrc)
outcomes <- list(
  OS  = list(end = "timeD",     event = "death_all"),
  PFS = list(end = "timeD_PFS", event = "PFS")
)

# Panel lettering is drawn INSIDE the top margin of each panel.
# This prevents A and B from being clipped off the top of the PDF.
add_panel_header <- function(letter, heading) {
  title(main = heading, line = 1.7, cex.main = 0.82, font.main = 2)
  mtext(letter, side = 3, line = 0.15, adj = 0,
        font = 2, cex = 1.10, col = "black")
}

format_ph_p <- function(p) {
  if (!is.finite(p)) return("P = NA")
  if (p < 0.0001) "P < 0.0001" else sprintf("P = %.4f", p)
}

plot_schoenfeld <- function(ph, group, outcome, panel_letter) {
  p <- ph$table["DNAmTL", "p"]
  label <- format_ph_p(p)

  # Use the actual fitted spline and bands returned by plot.cox.zph
  smooth <- plot(ph, var = "DNAmTL", se = TRUE, plot = FALSE)
  residuals <- ph$y[, "DNAmTL"]
  ok <- is.finite(ph$x) & is.finite(residuals)
  y_limits <- range(c(residuals[ok], smooth$y), finite = TRUE)
  extra <- diff(y_limits) * 0.04
  y_limits <- y_limits + c(-extra, extra)

  plot(
    ph, var = "DNAmTL", resid = FALSE, se = TRUE,
    col = c("red", "red"), lty = c(1, 2), lwd = c(2, 1),
    ylim = y_limits,
    xlab = "Follow-up time (days)",
    ylab = "Beta(t) for DNAmTL",
    main = ""
  )

  points(
    ph$x[ok], residuals[ok], pch = 1, cex = 0.50,
    col = adjustcolor("black", alpha.f = 0.50)
  )

  # Redraw the fitted curves on top of residual points
  lines(smooth$x, smooth$y[, 1], col = "red", lwd = 2)
  lines(smooth$x, smooth$y[, 2], col = "red", lwd = 1, lty = 2)
  lines(smooth$x, smooth$y[, 3], col = "red", lwd = 1, lty = 2)

  legend("topright", legend = label, bty = "n",
         cex = 0.78, inset = c(0.01, 0.01))

  add_panel_header(panel_letter,
                   paste0("Schoenfeld residuals (", group, ", ", outcome, ")"))
}

plot_loglog <- function(km_fit, group, outcome, panel_letter) {
  plot(
    km_fit, fun = "cloglog", col = c("red", "blue"), lwd = 1.7,
    xlab = "Follow-up time (days)",
    ylab = "Log[-log(Survival probability)]", main = ""
  )

  legend(
    "topleft", legend = c("Low DNAmTL", "High DNAmTL"),
    col = c("red", "blue"), lwd = 1.7,
    bty = "n", cex = 0.70, inset = c(0.01, 0.01)
  )

  add_panel_header(panel_letter,
                   paste0("Log-minus-log survival (", group, ", ", outcome, ")"))
}

# Open an A4 portrait, vector PDF
pdf(
  result_file("DNAmTL_PH_diagnostics_8panel_A4_final.pdf"),
  width = 8.27, height = 11.69, onefile = TRUE,
  family = "Helvetica", pointsize = 13
)

op <- par(no.readonly = TRUE)
par(
  mfrow = c(4, 2),
  mar = c(3.5, 4.0, 3.1, 0.5),
  oma = c(0.4, 0.3, 0.4, 0.3),
  mgp = c(2.1, 0.55, 0),
  cex = 0.90,
  cex.lab = 0.92,
  cex.axis = 0.85,
  tcl = -0.25
)

ph_tables <- list()
models <- list()
panel_number <- 0L

for (group in names(data_groups)) {
  for (outcome in names(outcomes)) {
    end_var <- outcomes[[outcome]]$end
    event_var <- outcomes[[outcome]]$event
    d <- data_groups[[group]]

    needed <- unique(c("DNAmTL", "time_blood", end_var,
                       event_var, adjust_vars))
    d <- d[complete.cases(d[, needed, drop = FALSE]), , drop = FALSE]

    # Valid delayed-entry intervals
    d <- d[is.finite(d$time_blood) & is.finite(d[[end_var]]) &
             d$time_blood >= 0 & d[[end_var]] > d$time_blood,
           , drop = FALSE]

    # Median grouping for graphical illustration only
    med <- median(d$DNAmTL, na.rm = TRUE)
    d$DNAmTL_dich <- factor(
      ifelse(d$DNAmTL > med, "high", "low"),
      levels = c("low", "high")
    )

    survival_formula <- as.formula(paste0(
      "Surv(time_blood, ", end_var, ", ", event_var, ") ~ ",
      paste(c("DNAmTL", adjust_vars), collapse = " + ")
    ))
    km_formula <- as.formula(paste0(
      "Surv(time_blood, ", end_var, ", ", event_var,
      ") ~ DNAmTL_dich"
    ))

    cox_model <- coxph(survival_formula, data = d, x = TRUE)
    ph <- cox.zph(cox_model, transform = "km", terms = TRUE, global = TRUE)
    km_fit <- survfit(km_formula, data = d)

    id <- paste(group, outcome, sep = "_")
    models[[id]] <- cox_model
    ph_tables[[id]] <- data.frame(
      Group = group, Outcome = outcome,
      Term = rownames(ph$table),
      Chi_square = ph$table[, "chisq"],
      DF = ph$table[, "df"],
      P_value = ph$table[, "p"],
      row.names = NULL
    )

    panel_number <- panel_number + 1L
    plot_schoenfeld(ph, group, outcome, LETTERS[panel_number])

    panel_number <- panel_number + 1L
    plot_loglog(km_fit, group, outcome, LETTERS[panel_number])
  }
}

par(op)
dev.off()

# Export all covariate-specific and global tests
ph_results <- do.call(rbind, ph_tables)
rownames(ph_results) <- NULL
print(ph_results)

write_results(ph_results, "PH_diagnostic_tests.xlsx")

# Convenient summary of the DNAmTL-specific PH tests
print(subset(ph_results, Term == "DNAmTL"))


# Extended Cox sensitivity analysis for LOCRC: DNAmTL x log(1 + time in years).
# The original model specification is retained; evidence is recomputed each run.
library(survival)

# Follow-up times (years)
years <- c(1, 3, 5, 10)

# 1. Define adjusted Cox model formula
covariates <- paste(
  c(
    "DNAmTL", "Age_diag", "Sex",
    "Schooling_yearsdich", "TNMdich",
    "locationdich", "chemradther",
    "cardiovad", "diabetesT2",
    "smokingdich", "lifetimealcohol_day",
    "tt(DNAmTL)"
  ),
  collapse = " + "
)

# 2. Function to fit extended Cox model
fit_extended_cox <- function(data, time, event) {

  formula <- as.formula(
    paste0(
      "Surv(time_blood, ", time, ", ", event,
      ") ~ ", covariates
    )
  )

  outcome <- if (event == "death_all") "OS" else "PFS"
  data <- model_data(data, formula, outcome)
  coxph(
    formula,
    data = data,
    tt = function(x, t, ...) {
      x * log1p(t / 365.25)
    }
  )
}

# Fit OS and PFS models
fit_os <- fit_extended_cox(
  locrc, "timeD", "death_all"
)

fit_pfs <- fit_extended_cox(
  locrc, "timeD_PFS", "PFS"
)


# 3. Function to calculate time-specific HRs (95% CI)
time_specific_hr <- function(model, years, outcome) {

  b <- coef(model)
  V <- vcov(model)

  # Coefficients for DNAmTL and its time interaction
  idx <- c("DNAmTL", "tt(DNAmTL)")

  if (!all(idx %in% names(b))) {
    stop("DNAmTL coefficients not found in the model.")
  }

  results <- lapply(years, function(year) {

    # Time transformation, consistent with model
    g <- log1p(year)

    L <- c(1, g)

    beta <- sum(L * b[idx])

    se <- sqrt(
      as.numeric(
        t(L) %*%
          V[idx, idx, drop = FALSE] %*%
          L
      )
    )

    HR <- exp(beta)
    CI_low <- exp(beta - 1.96 * se)
    CI_high <- exp(beta + 1.96 * se)

    data.frame(
      Outcome = outcome,
      Years = year,
      `HR (95% CI)` = sprintf(
        "%.2f (%.2f–%.2f)",
        HR, CI_low, CI_high
      ),
      check.names = FALSE
    )
  })

  do.call(rbind, results)
}


# 4. Calculate OS and PFS results
os_results <- time_specific_hr(
  fit_os, years, "Overall survival"
)

pfs_results <- time_specific_hr(
  fit_pfs, years, "Progression-free survival"
)


# 5. Combine into publication-ready dataframe
final_results <- rbind(
  os_results,
  pfs_results
)

rownames(final_results) <- NULL

print(final_results)


library(writexl)
write_xlsx(final_results,
           result_file("LOCRC_Time_Specific_HR.xlsx"))

save_session_info("04_PH_diagnostics")
