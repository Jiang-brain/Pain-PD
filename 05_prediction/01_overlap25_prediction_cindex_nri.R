source("00_setup/ukb_helpers.R")
source("00_setup/model_helpers.R")
load_packages(c("data.table", "survival"))

out_path <- result_path("overlap25_prediction_cindex_5y_10y_full_followup.csv")
nri_out_path <- result_path("overlap25_prediction_nri_5y_10y_full_followup.csv")

same <- data.table::fread(INTERSECTION_FILE)
same[, same_direction := sign(beta_pain_score) == sign(beta_PD_cox)]
overlap_proteins <- same[same_direction == TRUE, protein]
if (length(overlap_proteins) != 25) stop("Expected 25 same-direction proteins; found ", length(overlap_proteins), ".")

outcome_data <- make_outcome_data(UKB_DATA_FILE)
covars <- make_covariates(UKB_DATA_FILE, n_pcs = 10)
protein_cov <- data.table::fread(PROTEIN_COVARIATE_FILE)
data.table::setnames(protein_cov, "X74.0.0_zscore", "protein_cov_zscore")
protein_cov[, season_final := factor(season_final)]
protein_cov[, Batch := factor(Batch)]
protein_times <- data.table::fread(PROTEIN_TIME_FILE)
time_cols <- setdiff(names(protein_times), "eid")
protein_times[, protein_time := rowMeans(.SD, na.rm = TRUE), .SDcols = time_cols]
protein_times <- protein_times[, .(eid, protein_time)]

protein_raw <- data.table::fread(PROTEIN_MATRIX_FILE, select = c("eid", "ins_index", overlap_proteins))
protein_raw <- protein_raw[order(eid, ins_index)]
protein_raw <- protein_raw[!duplicated(eid)]
protein_raw[, ins_index := NULL]
protein_term_names <- make.unique(make.names(paste0("protein_", overlap_proteins)))
for (i in seq_along(overlap_proteins)) {
  protein <- overlap_proteins[i]
  values <- protein_raw[[protein]]
  protein_raw[[protein_term_names[i]]] <- as.numeric(scale(values))
}
protein_raw <- protein_raw[, c("eid", protein_term_names), with = FALSE]

analysis_data <- merge(outcome_data, covars, by = "eid")
analysis_data <- apply_main_exclusions(analysis_data)
analysis_data <- merge(as.data.table(analysis_data), protein_cov, by = "eid")
analysis_data <- merge(analysis_data, protein_times, by = "eid")
analysis_data <- merge(analysis_data, protein_raw, by = "eid")

covariates <- c(MAIN_COVARIATES, PROTEOMIC_TECHNICAL_COVARIATES)

fit_task <- function(task_name, task_data) {
  model_cols <- c("time_horizon", "status_horizon", covariates, protein_term_names)
  model_data <- as.data.frame(task_data[, model_cols, with = FALSE])
  model_data <- model_data[complete.cases(model_data), ]
  model_data <- droplevels(model_data)
  covs <- usable_covariates(model_data, covariates)
  cov_terms <- paste(covs, collapse = " + ")
  protein_terms <- paste(protein_term_names, collapse = " + ")

  base_fit <- survival::coxph(as.formula(paste("Surv(time_horizon, status_horizon) ~", cov_terms)), data = model_data)
  protein_alone_fit <- survival::coxph(as.formula(paste("Surv(time_horizon, status_horizon) ~", protein_terms)), data = model_data)
  protein_plus_fit <- survival::coxph(as.formula(paste("Surv(time_horizon, status_horizon) ~", cov_terms, "+", protein_terms)), data = model_data)

  cindex_rows <- data.table::rbindlist(list(
    cbind(task = task_name, cox_fit_row("traditional_covariates", base_fit, NULL, "traditional_covariates")),
    cbind(task = task_name, cox_fit_row("overlap25_proteins_alone", protein_alone_fit, NULL, "25_direction_aligned_overlap_proteins")),
    cbind(task = task_name, cox_fit_row("traditional_plus_overlap25_proteins", protein_plus_fit, base_fit, "25_direction_aligned_overlap_proteins"))
  ), fill = TRUE)
  eval_time <- max(model_data$time_horizon, na.rm = TRUE)
  nri <- continuous_nri_row(
    task_name,
    "traditional_covariates_vs_traditional_plus_overlap25_proteins",
    base_fit,
    protein_plus_fit,
    model_data,
    eval_time,
    status_col = "status_horizon"
  )
  list(cindex = cindex_rows, nri = nri)
}

task_specs <- list(
  list(name = "incident_PD_within_5y", data = make_horizon_task(analysis_data, 5)),
  list(name = "incident_PD_within_10y", data = make_horizon_task(analysis_data, 10)),
  list(name = "incident_PD_full_followup_after_2y", data = make_full_followup_task(analysis_data))
)

task_results <- lapply(task_specs, function(x) fit_task(x$name, x$data))
cindex_result <- data.table::rbindlist(lapply(task_results, `[[`, "cindex"), fill = TRUE)
nri_result <- data.table::rbindlist(lapply(task_results, `[[`, "nri"), fill = TRUE)

data.table::fwrite(cindex_result, out_path)
data.table::fwrite(nri_result, nri_out_path)
print(cindex_result)
print(nri_result)
cat("C-index output:", out_path, "\n")
cat("NRI output:", nri_out_path, "\n")

