source("00_setup/ukb_helpers.R")
load_packages(c("data.table"))

outcome <- make_outcome_data(UKB_DATA_FILE)
data.table::fwrite(outcome, result_path("ukb_followup_outcome.csv"))

cat("Rows:", nrow(outcome), "\n")
cat("PD events:", sum(outcome$status == 1, na.rm = TRUE), "\n")
cat("Output:", result_path("ukb_followup_outcome.csv"), "\n")

