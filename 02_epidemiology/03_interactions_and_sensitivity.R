source("00_setup/ukb_helpers.R")
load_packages(c("data.table", "survival"))

analysis_path <- result_path("ukb_main_analysis_dataset.csv")
if (!file.exists(analysis_path)) stop("Run 01_data_preparation/04_prepare_ukb_analysis_dataset.R first.")

dat <- data.table::fread(analysis_path)
effect_modifiers <- c("age", "gender", "ethnicity", "deprivation", "income", "university_education",
                      "smoking", "alcohol", "sedentary", "obesity", "prs")
main_covariates <- MAIN_COVARIATES

interaction_p <- function(exposure, modifier) {
  covariates <- setdiff(main_covariates, c(exposure, modifier))
  f0 <- as.formula(paste("Surv(time, status) ~", paste(c(exposure, modifier, covariates), collapse = " + ")))
  f1 <- as.formula(paste("Surv(time, status) ~", paste(c(paste0(exposure, "*", modifier), covariates), collapse = " + ")))
  fit0 <- survival::coxph(f0, data = dat)
  fit1 <- survival::coxph(f1, data = dat)
  anova(fit0, fit1, test = "Chisq")[2, "P(>|Chi|)"]
}

exposures <- c("acute_pain_status", "chronic_pain_status")
interaction_results <- expand.grid(exposure = exposures, modifier = effect_modifiers, stringsAsFactors = FALSE)
interaction_results$p_value <- mapply(interaction_p, interaction_results$exposure, interaction_results$modifier)
interaction_results$bonferroni_threshold <- 0.05 / nrow(interaction_results)
data.table::fwrite(interaction_results, result_path("ukb_pain_status_interactions.csv"))

# Extended-covariate sensitivity is listed here because the required source variables
# are not all available in the cleaned main dataset by default. Rebuild them upstream
# before enabling these models.
cat("Interaction output:", result_path("ukb_pain_status_interactions.csv"), "\n")
cat("Sensitivity analyses with dementia exclusion, pain medication, vascular disease, cancer, diabetes, and sleep require those variables in the input dataset.\n")

