source("00_setup/ukb_helpers.R")
load_packages(c("data.table"))

outcome_path <- result_path("ukb_followup_outcome.csv")
pain_path <- result_path("ukb_pain_exposures.csv")
covariate_path <- result_path("ukb_covariates.csv")

if (!file.exists(outcome_path)) stop("Run 01_data_preparation/01_build_ukb_followup_outcome.R first.")
if (!file.exists(pain_path)) stop("Run 01_data_preparation/02_define_ukb_pain_exposures.R first.")
if (!file.exists(covariate_path)) stop("Run 01_data_preparation/03_define_ukb_covariates.R first.")

outcome <- data.table::fread(outcome_path)
pain <- data.table::fread(pain_path)
covariates <- data.table::fread(covariate_path)

analysis_data <- merge(outcome, pain, by = "eid")
analysis_data <- merge(analysis_data, covariates, by = "eid")
analysis_data <- apply_main_exclusions(as.data.frame(analysis_data))

data.table::fwrite(as.data.table(analysis_data), result_path("ukb_main_analysis_dataset.csv"))

cat("Rows after complete-case, age >= 40, and 2-year incident PD exclusion:", nrow(analysis_data), "\n")
cat("PD events:", sum(analysis_data$status == 1, na.rm = TRUE), "\n")
cat("Output:", result_path("ukb_main_analysis_dataset.csv"), "\n")

