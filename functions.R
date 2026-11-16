# Shared utilities for the DNAmTL colorectal cancer analyses.
# Run scripts from the project root, for example source("5_cox_models.R").
# This file defines functions only: sourcing it does not load participant data.

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         ". See install_dependencies.R.", call. = FALSE)
  }
}

assert_columns <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(data)
}

# Default inputs are project-relative. An interactive file picker is offered
# when the analysis dataset has not been placed in data/. No personal path is
# stored in the scripts. An option supports batch execution with external data.
input_file <- function(name) {
  candidates <- c(file.path("data", name), name)
  found <- candidates[file.exists(candidates)]
  if (length(found)) return(found[1])
  if (interactive()) {
    message("Select ", name)
    return(file.choose())
  }
  stop("Input not found: data/", name, ". See data/README.md.", call. = FALSE)
}

result_file <- function(...) {
  path <- file.path(getOption("dnamtl.results_dir", "results"), ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

read_analysis_data <- function() {
  file <- getOption("dnamtl.data_file")
  if (is.null(file)) file <- input_file("df_followupdate.RData")
  # Load into a private environment, rather than replacing global R objects.
  env <- new.env(parent = emptyenv())
  load(file, envir = env)
  if (!exists("df_all", envir = env, inherits = FALSE)) {
    stop("The RData file must contain a data frame named df_all.")
  }
  data <- env$df_all
  stopifnot(is.data.frame(data))
  assert_columns(data, c("Age_diag", "DNAmTL"))
  if (!is.numeric(data$Age_diag) || !is.numeric(data$DNAmTL)) {
    stop("Age_diag and DNAmTL must be numeric.")
  }
  # Derive groups from age, avoiding fragile Unicode/string matching.
  data$age_group55 <- age_group(data$Age_diag, 55)
  data$age_group50 <- age_group(data$Age_diag, 50)
  data$newagesubgroup2 <- factor(
    data$age_group55, levels = c("EOCRC", "LOCRC"),
    labels = c("EOCRC (<55)", "LOCRC (>=55)")
  )
  # Retain existing coding; derive education only if absent.
  if (!"Schooling_yearsdich" %in% names(data)) {
    assert_columns(data, "Schooling_years")
    education <- as.character(data$Schooling_years)
    data$Schooling_yearsdich <- factor(
      ifelse(education == "<9", "<9 years",
             ifelse(education %in% c("9-10", ">10"), ">=9 years", NA)),
      levels = c("<9 years", ">=9 years")
    )
  }
  data
}

age_group <- function(age, cutoff = 55) {
  factor(ifelse(age < cutoff, "EOCRC", "LOCRC"), levels = c("EOCRC", "LOCRC"))
}

split_age_groups <- function(data, cutoff = 55) {
  split(data, age_group(data$Age_diag, cutoff), drop = FALSE)
}

# Primary adjustment set, retained from the supplied Cox functions.
cox_covariates <- c(
  "Age_diag", "Sex", "TNMdich", "locationdich", "chemradther",
  "cardiovad", "diabetesT2", "Schooling_yearsdich",
  "smokingdich", "lifetimealcohol_day"
)

# The original scripts apply this common restriction BEFORE within-group SDs
# and medians are calculated. Outcome-specific missingness is handled at fitting.
cohort_covariates <- setdiff(cox_covariates, "Age_diag")

complete_data <- function(data, variables) {
  assert_columns(data, variables)
  data[complete.cases(data[, variables, drop = FALSE]), , drop = FALSE]
}

standardize_scores <- function(data, vars = "DNAmTL") {
  assert_columns(data, vars)
  for (v in vars) {
    x <- data[[v]]
    if (!is.numeric(x)) stop(v, " must be numeric.")
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) stop("Cannot standardize ", v, ": SD is zero or undefined.")
    data[[v]] <- as.numeric((x - mean(x, na.rm = TRUE)) / s)
  }
  data
}

dichotomize_scores <- function(data, vars = "DNAmTL") {
  assert_columns(data, vars)
  for (v in vars) {
    median_value <- median(data[[v]], na.rm = TRUE)
    if (!is.finite(median_value)) stop("No finite median for ", v)
    data[[paste0(v, "_dich")]] <- factor(
      ifelse(data[[v]] > median_value, "high", "low"), levels = c("low", "high")
    )
  }
  data
}

prepare_cox_groups <- function(data, cutoff = 55, scores = "DNAmTL") {
  data <- complete_data(data, cohort_covariates)
  groups <- split_age_groups(data, cutoff)
  lapply(groups, function(d) standardize_scores(dichotomize_scores(d, scores), scores))
}

outcome_spec <- function(outcome) {
  specs <- list(
    OS = list(end = "timeD", event = "death_all", expression = "death_all"),
    PFS = list(end = "timeD_PFS", event = "PFS", expression = "PFS"),
    CSS = list(end = "timeD", event = "death_crccp", expression = "death_crccp == 1"),
    TTR = list(end = "timeD_recurr", event = "recurr_cp", expression = "recurr_cp == 1")
  )
  if (!outcome %in% names(specs)) stop("Outcome must be OS, PFS, CSS or TTR.")
  specs[[outcome]]
}

survival_formula <- function(outcome, terms, entrytime = "time_blood") {
  spec <- outcome_spec(outcome)
  as.formula(paste0(
    "survival::Surv(", entrytime, ", ", spec$end, ", ", spec$expression,
    ") ~ ", paste(terms, collapse = " + ")
  ))
}

# Exclude incomplete records and invalid delayed-entry intervals explicitly.
# Never convert competing events to a multi-state response: CSS/TTR use == 1.
model_data <- function(data, formula, outcome, entrytime = "time_blood") {
  d <- complete_data(data, all.vars(formula))
  spec <- outcome_spec(outcome)
  valid <- is.finite(d[[entrytime]]) & is.finite(d[[spec$end]]) &
    d[[entrytime]] >= 0 & d[[spec$end]] > d[[entrytime]]
  if (any(!valid)) warning(sum(!valid), " invalid delayed-entry intervals excluded for ", outcome)
  d <- d[valid, , drop = FALSE]
  if (outcome %in% c("OS", "PFS") && !all(d[[spec$event]] %in% c(0, 1))) {
    stop("OS and PFS event variables must be coded 0/1.")
  }
  if (!nrow(d) || !any(d[[spec$event]] == 1)) stop("No usable events for ", outcome)
  d
}

format_p <- function(p, digits = 3, threshold = 0.001) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < threshold, paste0("<", threshold), formatC(p, digits = digits, format = "f")))
}

fit_cox_model <- function(data, score, outcome, entrytime = "time_blood",
                          covariates = cox_covariates, extra_covariates = character()) {
  require_packages("survival")
  adjustment <- unique(c(covariates, extra_covariates))
  formula <- survival_formula(outcome, c(score, adjustment), entrytime)
  d <- model_data(data, formula, outcome, entrytime)
  # A subgroup-defining variable (e.g. sex) cannot also be estimated within
  # that subgroup. Record constant adjustment variables rather than failing.
  constant <- adjustment[vapply(d[, adjustment, drop = FALSE],
                               function(x) length(unique(x)) < 2, logical(1))]
  adjustment <- setdiff(adjustment, constant)
  formula <- survival_formula(outcome, c(score, adjustment), entrytime)
  fit_warnings <- character()
  fit <- withCallingHandlers(
    survival::coxph(formula, data = d, x = TRUE, y = TRUE, model = TRUE),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  attr(fit, "excluded_constant_covariates") <- constant
  attr(fit, "fit_warnings") <- unique(fit_warnings)
  fit
}

cox_result <- function(fit, score, outcome, n_input) {
  tab <- summary(fit)$coefficients
  # Use model-matrix term assignments, not a substring grep that could match
  # a different score with a shared prefix.
  idx <- fit$assign[[score]]
  if (!length(idx)) stop("Score term not found in fitted model: ", score)
  term <- names(coef(fit))[idx]
  ci <- exp(confint(fit)[term, , drop = FALSE])
  hr <- exp(coef(fit)[term])
  p <- tab[term, "Pr(>|z|)"]
  data.frame(
    score = term, exposure = score, outcome = outcome,
    N_input = n_input, N = fit$n, N_events = fit$nevent,
    HR = unname(hr), LL = ci[, 1], UL = ci[, 2], P_value = unname(p),
    `aHR(95% CI)` = sprintf("%.2f (%.2f, %.2f)", hr, ci[, 1], ci[, 2]),
    aPvalue = format_p(p),
    omitted_covariates = paste(attr(fit, "excluded_constant_covariates"), collapse = "; "),
    Fit_status = if (length(attr(fit, "fit_warnings"))) "REVIEW_WARNINGS" else "OK",
    Fit_warnings = paste(attr(fit, "fit_warnings"), collapse = "; "),
    row.names = NULL, check.names = FALSE
  )
}

multicoxfull <- function(data, entrytime = "time_blood", score, outcome,
                         extra_covariates = character()) {
  fit <- fit_cox_model(data, score, outcome, entrytime, extra_covariates = extra_covariates)
  cox_result(fit, score, outcome, nrow(data))
}

multicoxfulltimesurg <- function(data, entrytime = "time_blood", score, outcome) {
  multicoxfull(data, entrytime, score, outcome, extra_covariates = "brelsurgtim_month")
}

# Present in the original external cox_complete.R, but missing from functions.R.
multicoxfulltimechemo <- function(data, entrytime = "time_blood", score, outcome) {
  multicoxfull(data, entrytime, score, outcome, extra_covariates = "anychemorebldt")
}

run_cox_grid <- function(groups, scores = c("DNAmTL", "DNAmTL_dich"),
                          outcomes = c("OS", "CSS", "PFS", "TTR"),
                          extra_covariates = character()) {
  results <- list()
  for (group in names(groups)) {
    for (outcome in outcomes) {
      for (score in scores) {
        result <- multicoxfull(groups[[group]], score = score, outcome = outcome,
                               extra_covariates = extra_covariates)
        result$dataset <- group
        results[[paste(group, outcome, score, sep = "_")]] <- result
      }
    }
  }
  result <- do.call(rbind, results)
  rownames(result) <- NULL
  n_flagged <- sum(result$Fit_status != "OK")
  if (n_flagged) warning(n_flagged, " fitted result rows need review; see Fit_warnings in the export.")
  result
}

write_results <- function(data, filename) {
  require_packages("writexl")
  writexl::write_xlsx(data, result_file(filename))
}

# Robust regression: coefficient CIs use the original t approximation;
# global factor P values use car::Anova's Type II Wald test.
multivariable_rlm_single <- function(data, y_vars, x_vars) {
  require_packages(c("MASS", "car"))
  results <- list()
  for (y in names(y_vars)) {
    d <- complete_data(data, c(y, x_vars))
    fit <- MASS::rlm(reformulate(x_vars, y), data = d, maxit = 100)
    if (!isTRUE(fit$converged)) stop("Robust regression did not converge for ", y)
    s <- summary(fit)
    global <- as.data.frame(car::Anova(fit, type = "II"))
    mm <- model.matrix(fit)
    assignment <- attr(mm, "assign")
    term_labels <- attr(terms(fit), "term.labels")
    for (term in setdiff(rownames(s$coefficients), "(Intercept)")) {
      base <- term_labels[assignment[match(term, colnames(mm))]]
      estimate <- s$coefficients[term, "Value"]
      se <- s$coefficients[term, "Std. Error"]
      critical <- qt(0.975, s$df[2])
      ci <- estimate + c(-1, 1) * critical * se
      results[[paste(y, term)]] <- data.frame(
        Y_Variable = y, Y_Label = y_vars[[y]], Term = term, Base_Variable = base,
        N = nrow(d), Estimate = estimate, Std_Error = se, CI_lower = ci[1], CI_upper = ci[2],
        estimateCI = sprintf("%.3f (%.3f, %.3f)", estimate, ci[1], ci[2]),
        p_value = 2 * pt(-abs(estimate / se), s$df[2]),
        global_p_value = global[base, ncol(global)], row.names = NULL
      )
    }
  }
  do.call(rbind, results)
}

interaction_lrt <- function(data, modifier, outcome, score = "DNAmTL",
                            covariates = cox_covariates) {
  require_packages("survival")
  adjustment <- setdiff(covariates, c(score, modifier))
  f0 <- survival_formula(outcome, c(score, modifier, adjustment))
  f1 <- survival_formula(outcome, c(paste(score, modifier, sep = " * "), adjustment))
  # Both nested models must use exactly the same records.
  d <- model_data(data, f1, outcome)
  fit0 <- survival::coxph(f0, data = d)
  fit1 <- survival::coxph(f1, data = d)
  statistic <- 2 * (as.numeric(logLik(fit1)) - as.numeric(logLik(fit0)))
  degrees <- attr(logLik(fit1), "df") - attr(logLik(fit0), "df")
  data.frame(outcome = outcome, exposure = score, modifier = modifier, N = fit1$n,
             Chi_square = statistic, DF = degrees,
             P_value = pchisq(statistic, degrees, lower.tail = FALSE))
}

get_cindex_results <- function(data, time_var, event_var, group_name, outcome_name) {
  require_packages(c("survival", "compareC"))
  variables <- c(time_var, event_var, "Age_diag", "Sex", "TNM_adj", "DNAmTL")
  d <- complete_data(data, variables)
  if (!all(d[[event_var]] %in% c(0, 1))) stop("C-index event must be coded 0/1.")
  response <- paste0("survival::Surv(", time_var, ", ", event_var, ")")
  fit_base <- survival::coxph(as.formula(paste(response, "~ Age_diag + Sex + TNM_adj")), data = d)
  fit_mtl <- survival::coxph(as.formula(paste(response, "~ Age_diag + Sex + TNM_adj + DNAmTL")), data = d)
  base <- survival::concordance(fit_base)
  mtl <- survival::concordance(fit_mtl)
  # compareC expects scores increasing with survival, whereas Cox linear
  # predictors increase with hazard. Negate both predictors to align direction.
  comparison <- compareC::compareC(d[[time_var]], d[[event_var]],
                                  scoreY = -predict(fit_base), scoreZ = -predict(fit_mtl))
  format_c <- function(x) {
    ci <- x$concordance + c(-1, 1) * qnorm(0.975) * sqrt(x$var)
    sprintf("%.2f (%.2f-%.2f)", x$concordance, ci[1], ci[2])
  }
  data.frame(
    Group = group_name, Outcome = outcome_name, N = nrow(d),
    Model1_Base_Cindex = format_c(base), Model2_Add_DNAmTL_Cindex = format_c(mtl),
    C_base = base$concordance, C_DNAmTL = mtl$concordance,
    Delta_C = mtl$concordance - base$concordance,
    P_value = comparison$pval, P_value_Comparison = format_p(comparison$pval, 2)
  )
}

plot_DNAatri_comparisonnew <- function(data, score_var = "TL_shortenmeanbirh",
                                      group_var, xtitle) {
  require_packages(c("ggplot2", "ggpubr"))
  d <- complete_data(data, c(score_var, group_var, "newagesubgroup2"))
  ggplot2::ggplot(d, ggplot2::aes(x = .data[[group_var]], y = .data[[score_var]],
                                 fill = newagesubgroup2)) +
    ggplot2::geom_boxplot(width = 0.6, outlier.shape = NA) +
    ggplot2::scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
    ggplot2::labs(x = xtitle, y = "DNAmTL shortening index", fill = NULL) +
    ggpubr::stat_compare_means(ggplot2::aes(group = newagesubgroup2),
                              method = "kruskal.test", label = "p.signif") +
    ggplot2::theme_minimal() + ggplot2::theme(legend.position = "top")
}

save_session_info <- function(label) {
  writeLines(capture.output(sessionInfo()), result_file(paste0(label, "_sessionInfo.txt")))
}
