# Baseline tables, blood-collection timing, and descriptive DNAmTL figures.
# Input: df_all in df_followupdate.RData. Outputs are written to results/.
source("functions.R", local = TRUE)
require_packages(c("tableone", "ggplot2", "ggpubr", "ggpmisc"))
library(ggplot2)
df_all <- read_analysis_data()

baseline_vars <- c(
  "Age_diag", "Sex", "Schooling_years", "smoking", "liftimdrink_cat",
  "BMI", "BMI_cat", "BMI_5_14earlier", "BMI5_14er_cat", "cardiovad",
  "diabetesT2", "high_blood_pressure", "other_cancer", "TNM_adj", "Location",
  "chemradther", "neotreat", "MSI_gent", "brafmut", "krasmut", "cimphi"
)
blood_vars <- c("neotreat", "chemradther", "brelsurgtim_month",
                "brelsurgtim_monthcat", "anychemorebldt", "anychemocat")

make_table <- function(data, variables, nonnormal, include_missing, test = TRUE) {
  assert_columns(data, variables)
  tab <- tableone::CreateTableOne(vars = variables, strata = "newagesubgroup2",
                                  data = data, includeNA = include_missing, test = test)
  printed <- print(tab, nonnormal = nonnormal, catDigits = 1, contDigits = 1,
                   showAllLevels = TRUE, missing = include_missing, quote = FALSE,
                   noSpaces = TRUE, printToggle = FALSE)
  data.frame(Characteristic = trimws(rownames(printed)), printed,
             row.names = NULL, check.names = FALSE)
}
write_results(make_table(df_all, baseline_vars, c("Age_diag", "BMI", "BMI_5_14earlier"), FALSE),
              "table1nomissing.xlsx")
write_results(make_table(df_all, baseline_vars, c("Age_diag", "BMI", "BMI_5_14earlier"), TRUE, FALSE),
              "table1withmissing.xlsx")
write_results(list(
  without_missing = make_table(df_all, blood_vars, "brelsurgtim_month", FALSE),
  with_missing = make_table(df_all, blood_vars, "brelsurgtim_month", TRUE)
), "table1bloodtime.xlsx")

# Compute correlations from the current data instead of embedding old results.
groups <- c(list(All = df_all), split_age_groups(df_all))
correlations <- do.call(rbind, lapply(names(groups), function(group) {
  d <- complete_data(groups[[group]], c("Age_diag", "DNAmTL"))
  test <- cor.test(d$Age_diag, d$DNAmTL, method = "pearson")
  data.frame(Group = group, N = nrow(d), Pearson_r = unname(test$estimate), P_value = test$p.value)
}))
write_results(correlations, "DNAmTL_age_correlations.xlsx")

scatter <- ggplot(df_all, aes(Age_diag, DNAmTL, colour = newagesubgroup2)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  ggpmisc::stat_poly_eq(aes(label = after_stat(eq.label)), formula = y ~ x,
                       parse = TRUE, size = 3.5) +
  scale_colour_manual(values = c("#E41A1C", "#377EB8")) +
  labs(x = "Age at diagnosis (years)", y = "DNAmTL", colour = NULL) +
  theme_classic() + theme(legend.position = "top")

assert_columns(df_all, "TL_shortenmeanbirh")
# The supplied index is retained; it is descriptive, not a longitudinal measure.
boxplot <- ggplot(df_all, aes(newagesubgroup2, TL_shortenmeanbirh, fill = newagesubgroup2)) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  ggpubr::stat_compare_means(method = "kruskal.test", label = "p.signif") +
  labs(x = NULL, y = "DNAmTL shortening index") +
  theme_minimal() + theme(legend.position = "none")
figure <- ggpubr::ggarrange(scatter, boxplot, ncol = 2, labels = c("A", "B"))
ggsave(result_file("Fig1_scatter_box.pdf"), figure, width = 10, height = 5)

group_vars <- c(Sex = "Sex", smoking = "Smoking status",
                lifealcoholdich = "Lifetime alcohol consumption", TNM_adj = "TNM stage",
                Location = "Tumour location", chemradther = "Adjuvant treatment")
plots <- lapply(names(group_vars), function(v) {
  plot_DNAatri_comparisonnew(df_all, group_var = v, xtitle = group_vars[[v]])
})
figure <- ggpubr::ggarrange(plotlist = plots, ncol = 3, nrow = 2,
                           common.legend = TRUE, legend = "top")
ggsave(result_file("shortening_index_subgroups.pdf"), figure, width = 12, height = 8)

# Formal slope comparison: age-by-onset-group interaction.
slope_fit <- lm(DNAmTL ~ Age_diag * age_group55, data = df_all)
tab <- as.data.frame(coef(summary(slope_fit)))
tab$Term <- rownames(tab)
write_results(tab, "DNAmTL_age_slope_interaction.xlsx")
save_session_info("02_descriptive")
