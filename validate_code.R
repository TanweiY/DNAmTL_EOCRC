# Small reproducibility checks that do not require participant-level data.
# Run from the project root: Rscript validate_code.R
for (file in list.files(pattern = "\\.R$")) parse(file, encoding = "UTF-8")
source("functions.R", local = TRUE)

# Boundary cases: no participant aged 50-54 is lost at either cutoff.
ages <- c(49, 50, 54, 55, NA)
stopifnot(identical(as.character(age_group(ages, 50)), c("EOCRC", "LOCRC", "LOCRC", "LOCRC", NA_character_)))
stopifnot(identical(as.character(age_group(ages, 55)), c("EOCRC", "EOCRC", "EOCRC", "LOCRC", NA_character_)))

# Check that input loading is isolated and does not change the original object.
temporary <- tempfile(fileext = ".RData")
df_all <- data.frame(Age_diag = ages, DNAmTL = c(6.8, 6.7, 6.6, 6.5, NA),
                      Schooling_years = c("<9", "9-10", ">10", "<9", NA))
original <- df_all
save(df_all, file = temporary)
old <- options(dnamtl.data_file = temporary)
loaded <- read_analysis_data()
options(old)
stopifnot(identical(df_all, original), nrow(loaded) == nrow(original))
unlink(temporary)

# A model must report fitted N, excluding missing predictors and invalid intervals.
require_packages(c("survival", "compareC"))
set.seed(20261116)
n <- 300
data <- data.frame(Age_diag = runif(n, 30, 85), Sex = factor(sample(c("Female", "Male"), n, TRUE)),
                    TNM_adj = factor(sample(c("I", "II", "III", "IV"), n, TRUE)),
                    DNAmTL = rnorm(n), time_blood = rep(0, n))
event_time <- rexp(n, exp(-0.8 * data$DNAmTL) / 100)
censor_time <- rexp(n, 1 / 150)
data$timeD <- pmin(event_time, censor_time)
data$death_all <- as.integer(event_time <= censor_time)
data$DNAmTL[1] <- NA
data$time_blood[2] <- data$timeD[2]
fit <- suppressWarnings(fit_cox_model(data, "DNAmTL", "OS", covariates = c("Age_diag", "Sex", "TNM_adj")))
result <- cox_result(fit, "DNAmTL", "OS", n)
stopifnot(result$N_input == n, result$N == n - 2, result$HR < 1)

# The survival-score direction passed to compareC must agree with concordance().
paired <- complete_data(data, c("timeD", "death_all", "Age_diag", "Sex", "TNM_adj", "DNAmTL"))
base <- survival::coxph(survival::Surv(timeD, death_all) ~ Age_diag + Sex + TNM_adj, data = paired)
extended <- survival::coxph(survival::Surv(timeD, death_all) ~ Age_diag + Sex + TNM_adj + DNAmTL, data = paired)
comparison <- compareC::compareC(paired$timeD, paired$death_all, -predict(base), -predict(extended))
stopifnot(abs(unname(comparison$est.c[1]) - survival::concordance(base)$concordance) < 1e-8,
          abs(unname(comparison$est.c[2]) - survival::concordance(extended)$concordance) < 1e-8)

# A p-value rounded to zero must never be described as below a smaller threshold.
stopifnot(format_p(0.0004) == "<0.001", format_p(0.0124) == "0.012")
cat("All syntax and synthetic-data checks passed.\n")
