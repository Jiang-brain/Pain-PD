source("00_setup/ukb_helpers.R")
load_packages(c("data.table"))

covariates <- make_covariates(UKB_DATA_FILE, n_pcs = 10)
data.table::fwrite(as.data.table(covariates), result_path("ukb_covariates.csv"))

cat("Rows:", nrow(covariates), "\n")
cat("Output:", result_path("ukb_covariates.csv"), "\n")

